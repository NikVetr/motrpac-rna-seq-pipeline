#!/usr/bin/env python3
"""Capture task inputs, counters and resources without downloading large data files."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import subprocess

from summarize_workflow_cost import monitor


SMALL_METRICS = (".umi_metrics.json", ".umi_molecule_expression_metrics.json",
                 ".sampling_manifest.json", ".Log.final.out", "_qc_info.csv", ".cnt",
                 "expression_metadata.tsv")


def cloud(*args):
    try:
        return subprocess.check_output(["gcloud", *args], stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as error:
        raise ValueError(error.stderr.decode(errors="replace").strip()) from error


def file_uris(value):
    if isinstance(value, str) and value.startswith("gs://"):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from file_uris(child)
    elif isinstance(value, list):
        for child in value:
            yield from file_uris(child)


def collect(metadata_path, output, workers=8):
    metadata = json.loads(metadata_path.read_text())
    if not metadata.get("calls"):
        raise ValueError("metadata must include calls, inputs, outputs and runtime attributes")
    output.mkdir(parents=True, exist_ok=False)
    (output / "metadata.json").write_bytes(metadata_path.read_bytes())
    workflow_inputs = metadata.get("inputs") or json.loads(
        metadata.get("submittedFiles", {}).get("inputs", "{}"))
    attempts = [(call, attempt) for call, calls in metadata["calls"].items() for attempt in calls]
    # No running-stream snapshot is presented as a completed resource measurement.
    terminal = {"Done", "Failed", "RetryableFailure", "Aborted", "Bypassed"}
    selected = [(call, attempt) for call, attempt in attempts
                if attempt.get("executionStatus") in terminal]
    uris = sorted({uri for _, attempt in selected
                   for field in ("inputs", "outputs", "monitoringLog", "stdout", "stderr")
                   for uri in file_uris(attempt.get(field))})

    def describe(uri):
        obj = json.loads(cloud("storage", "objects", "describe", uri, "--format=json"))
        size = int(obj.get("size", obj.get("size_bytes", -1)))
        if size < 0 or not str(obj.get("generation", "")).isdigit():
            raise ValueError(f"object lacks size/generation: {uri}")
        return uri, {"size_bytes": size, "generation": str(obj["generation"])}

    objects, errors = {}, []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [(uri, pool.submit(describe, uri)) for uri in uris]
        for uri, future in futures:
            try:
                key, value = future.result()
                objects[key] = value
            except (ValueError, subprocess.CalledProcessError) as error:
                errors.append({"uri": uri, "error": str(error)})
    (output / "objects.json").write_text(json.dumps(objects, indent=2) + "\n")

    def download(uri, target, limit):
        obj = objects.get(uri)
        if obj is None:
            return None
        if obj["size_bytes"] > limit:
            raise ValueError(f"profiling artifact exceeds {limit} bytes: {uri}")
        versioned = uri.split("#", 1)[0] + "#" + obj["generation"]
        data = cloud("storage", "cat", versioned)
        if len(data) != obj["size_bytes"]:
            raise ValueError(f"profiling artifact size changed: {uri}")
        target.write_bytes(data)
        return data

    rows = []
    for call, attempt in selected:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", call):
            raise ValueError(f"invalid call name: {call}")
        shard, number = int(attempt["shardIndex"]), int(attempt["attempt"])
        directory = output / f"{call}.shard-{shard}.attempt-{number}"
        directory.mkdir()
        inputs = attempt.get("inputs", {})
        if attempt.get("jobId") and not inputs:
            raise ValueError(f"metadata lacks task inputs: {call}")
        samples = workflow_inputs.get("rnaseq_pipeline.sample_prefix", [])
        row = {"call": call, "shard": shard, "attempt": number,
               "start": attempt.get("start"), "end": attempt.get("end"),
               "sample": inputs.get("SID", inputs.get("sample_prefix",
                         samples[shard] if 0 <= shard < len(samples) else None)),
               "reference_release": workflow_inputs.get("rnaseq_pipeline.reference_release"),
               "status": attempt["executionStatus"], "inputs": inputs,
               "runtime": attempt.get("runtimeAttributes", {}),
               "cache": attempt.get("callCaching", {}), "failures": attempt.get("failures", []),
               "input_bytes": {key: sum(objects[uri]["size_bytes"] for uri in set(file_uris(value)))
                               for key, value in inputs.items()
                               if set(file_uris(value)) and all(uri in objects for uri in file_uris(value))},
               "monitoring": None}
        if attempt.get("jobId"):
            job = json.loads(cloud("batch", "jobs", "describe", attempt["jobId"], "--format=json"))
            (directory / "batch.json").write_text(json.dumps(job, indent=2) + "\n")
            row["batch"] = job
        for field in ("monitoringLog", "stdout", "stderr"):
            uri = attempt.get(field)
            if not uri:
                continue
            path = directory / field
            data = download(uri, path, 64 * 1024**2)
            if field == "monitoringLog" and data:
                row["monitoring"] = monitor(path)
        row["metrics"] = {}
        for uri in sorted(set(file_uris(attempt.get("outputs", {})))):
            if not uri.endswith(SMALL_METRICS):
                continue
            name = uri.rsplit("/", 1)[-1]
            target = directory / name
            if target.exists():
                raise ValueError(f"duplicate profiling output name: {name}")
            data = download(uri, target, 1024**2)
            if data is not None:
                row["metrics"][name] = json.loads(data) if name.endswith(".json") else data.decode()
        row["outputs"] = attempt.get("outputs", {})  # Includes scalar post-trim read counts.
        rows.append(row)
    result = {"schema_version": 1, "captured_utc": datetime.now(timezone.utc).isoformat(),
              "workflow_id": metadata["id"], "workflow_status": metadata.get("status"),
              "workflow_inputs": workflow_inputs,
              "submitted_files": metadata.get("submittedFiles", {}),
              "not_terminal": [{"call": c, "shard": a.get("shardIndex"), "attempt": a.get("attempt"),
                                "status": a.get("executionStatus")} for c, a in attempts
                               if a.get("executionStatus") not in terminal],
              "attempts": rows, "errors": errors, "complete_capture": not errors}
    def encode(value):
        if isinstance(value, Decimal):
            return float(value)
        raise TypeError(type(value).__name__)
    (output / "profiles.json").write_text(json.dumps(result, indent=2, default=encode) + "\n")
    manifest = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  ./{path.relative_to(output)}"
                for path in sorted(output.rglob("*")) if path.is_file()]
    (output / "evidence-manifest.sha256").write_text("\n".join(manifest) + "\n")
    if errors:
        raise ValueError(f"incomplete capture: {len(errors)} unavailable objects; see profiles.json")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("output", type=Path, help="new directory for this capture")
    args = parser.parse_args()
    result = collect(args.metadata, args.output)
    print(f"Captured {len(result['attempts'])} terminal attempts; {len(result['not_terminal'])} still pending")
