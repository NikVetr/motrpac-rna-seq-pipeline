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


RUN_ROOT = re.compile(r"^(gs://[^/]+/.+?/rnaseq_pipeline/[^/]+)/")
TERMINAL_WORKFLOW_STATES = {"Succeeded", "Failed", "Aborted"}
TERMINAL_ATTEMPT_STATES = {"Done", "Failed", "RetryableFailure", "Aborted", "Bypassed"}
SMALL_METRICS = (".umi_metrics.json", ".umi_molecule_expression_metrics.json",
                 ".sampling_manifest.json", ".Log.final.out", "_qc_info.csv", ".cnt",
                 "_contamination_sampling.json", ".qc_diagnostics.json", "expression_metadata.tsv",
                 ".rsem_convergence.tsv", ".rsem_gene_convergence.tsv")


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


def collect(metadata_path, output, workers=8, batch=True, cloud_read=None):
    if workers < 1:
        raise ValueError("workers must be positive")
    read = cloud_read or cloud
    metadata = json.loads(metadata_path.read_text())
    if not metadata.get("calls"):
        raise ValueError("metadata must include calls, inputs, outputs and runtime attributes")
    output.mkdir(parents=True, exist_ok=False)
    (output / "metadata.json").write_bytes(metadata_path.read_bytes())
    workflow_inputs = metadata.get("inputs") or json.loads(
        metadata.get("submittedFiles", {}).get("inputs", "{}"))
    attempts = [(call, attempt) for call, calls in metadata["calls"].items() for attempt in calls]
    # No running-stream snapshot is presented as a completed resource measurement.
    selected = [(call, attempt) for call, attempt in attempts
                if attempt.get("executionStatus") in TERMINAL_ATTEMPT_STATES]
    preempted_logs = {attempt.get("monitoringLog") for _, attempt in selected
                     if attempt.get("executionStatus") == "RetryableFailure"
                     and "VMPreemption(50001)" in json.dumps(attempt.get("failures", []))}
    preempted_logs -= {attempt.get("monitoringLog") for _, attempt in selected
                      if attempt.get("executionStatus") == "Done"}
    # Cache hits reuse another call's results without running, so they have no monitoring log, stdout or stderr
    # (those stay with the call that produced the results). Attempts that never finished (preempted, failed or
    # aborted) may also lack them. Only successful attempts that ran must have their logs.
    log_fields = ("monitoringLog", "stdout", "stderr")
    cached = lambda attempt: (attempt.get("callCaching") or {}).get("hit")
    cached_logs = {attempt.get(field) for _, attempt in selected if cached(attempt) for field in log_fields}
    unfinished_logs = {attempt.get(field) for _, attempt in selected
                       if attempt.get("executionStatus") != "Done" for field in log_fields}
    unfinished_logs -= {attempt.get(field) for _, attempt in selected
                        if attempt.get("executionStatus") == "Done" and not cached(attempt) for field in log_fields}
    uris = sorted({uri for _, attempt in selected
                   for field in ("inputs", "outputs", "monitoringLog", "stdout", "stderr")
                   for uri in file_uris(attempt.get(field))})

    def describe(uri):
        obj = json.loads(read("storage", "objects", "describe", uri, "--format=json"))
        size = int(obj.get("size", obj.get("size_bytes", -1)))
        if size < 0 or not str(obj.get("generation", "")).isdigit():
            raise ValueError(f"object lacks size/generation: {uri}")
        return uri, {"size_bytes": size, "generation": str(obj["generation"])}

    def list_root(root):
        """Size and generation of every object under one workflow root, in one paginated listing."""
        try:
            listing = json.loads(read("storage", "objects", "list", f"{root}/**", "--format=json(bucket,name,size,generation)"))
        except ValueError as error:
            if "matched no objects" not in str(error):
                raise
            listing = []
        return {f"gs://{obj['bucket']}/{obj['name']}": {"size_bytes": int(obj["size"]), "generation": str(obj["generation"])}
                for obj in listing}

    # Workflow files are listed per root (one call instead of one describe per object); others are described.
    roots = sorted({match[1] for match in map(RUN_ROOT.match, uris) if match})
    with ThreadPoolExecutor(max_workers=workers) as pool:
        listed = {key: value for part in pool.map(list_root, roots) for key, value in part.items()}

    def lookup(uri):
        match = RUN_ROOT.match(uri)
        if match is None:
            return describe(uri)
        if uri not in listed:
            raise ValueError(f"{uri} not found: 404 (absent from the listing of {match[1]})")
        return uri, listed[uri]

    objects, errors = {}, []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [(uri, pool.submit(lookup, uri)) for uri in uris]
        for uri, future in futures:
            try:
                key, value = future.result()
                objects[key] = value
            except (ValueError, subprocess.CalledProcessError) as error:
                missing = bool(re.search(r"not found:\s*404\b|\b404[:\s]+not found\b", str(error), re.I))
                errors.append({"uri": uri, "error": str(error),
                               "expected_after_preemption": uri in preempted_logs and missing,
                               "expected_after_cache_hit": uri in cached_logs and missing,
                               "expected_after_unfinished_attempt": uri in unfinished_logs and missing})
    (output / "objects.json").write_text(json.dumps(objects, indent=2) + "\n")

    def download(uri, target, limit):
        obj = objects.get(uri)
        if obj is None:
            return None
        if obj["size_bytes"] > limit:
            raise ValueError(f"profiling artifact exceeds {limit} bytes: {uri}")
        versioned = uri.split("#", 1)[0] + "#" + obj["generation"]
        data = read("storage", "cat", versioned)
        if len(data) != obj["size_bytes"]:
            raise ValueError(f"profiling artifact size changed: {uri}")
        target.write_bytes(data)
        return data

    def capture_attempt(call, attempt):
        """Write one attempt's directory; return its row and any Batch read errors."""
        attempt_errors = []
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", call):
            raise ValueError(f"invalid call name: {call}")
        shard, number = int(attempt["shardIndex"]), int(attempt["attempt"])
        directory = output / f"{call}.shard-{shard}.attempt-{number}"
        directory.mkdir()
        inputs = attempt.get("inputs", {})
        if attempt.get("jobId") and not inputs:
            raise ValueError(f"metadata lacks task inputs: {call}")
        samples = workflow_inputs.get("sample_prefix", workflow_inputs.get("rnaseq_pipeline.sample_prefix", []))
        row = {"call": call, "shard": shard, "attempt": number,
               "start": attempt.get("start"), "end": attempt.get("end"),
               "sample": inputs.get("SID", inputs.get("sample_prefix",
                         samples[shard] if 0 <= shard < len(samples) else None)),
               "reference_release": workflow_inputs.get("reference_release", workflow_inputs.get("rnaseq_pipeline.reference_release")),
               "status": attempt["executionStatus"], "inputs": inputs,
               "runtime": attempt.get("runtimeAttributes", {}),
               "cache": attempt.get("callCaching", {}), "failures": attempt.get("failures", []),
               "input_bytes": {key: sum(objects[uri]["size_bytes"] for uri in set(file_uris(value)))
                               for key, value in inputs.items()
                               if set(file_uris(value)) and all(uri in objects for uri in file_uris(value))},
               "monitoring": None}
        if batch and attempt.get("jobId"):
            try:
                job = json.loads(read("batch", "jobs", "describe", attempt["jobId"], "--format=json"))
                (directory / "batch.json").write_text(json.dumps(job, indent=2) + "\n")
                row["batch"] = job
            except ValueError as error:
                attempt_errors.append({"job": attempt["jobId"], "error": str(error)})
        for field in log_fields:
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
            limit = 64 * 1024**2 if name.endswith((".rsem_convergence.tsv", ".rsem_gene_convergence.tsv")) else 1024**2
            data = download(uri, target, limit)
            if data is not None:
                row["metrics"][name] = json.loads(data) if name.endswith(".json") else data.decode()
        row["outputs"] = attempt.get("outputs", {})  # Includes scalar post-trim read counts.
        return row, attempt_errors

    # Bound concurrent network reads independently of the number of task attempts.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        captured = list(pool.map(lambda pair: capture_attempt(*pair), selected))
    rows = [row for row, _ in captured]
    errors += [error for _, attempt_errors in captured for error in attempt_errors]
    pending = [{"call": c, "shard": a.get("shardIndex"), "attempt": a.get("attempt"),
                "status": a.get("executionStatus")} for c, a in attempts
               if a.get("executionStatus") not in TERMINAL_ATTEMPT_STATES]
    result = {"schema_version": 1, "captured_utc": datetime.now(timezone.utc).isoformat(),
              "workflow_id": metadata["id"], "workflow_status": metadata.get("status"),
              "workflow_inputs": workflow_inputs,
              "submitted_files": metadata.get("submittedFiles", {}),
              "not_terminal": pending,
              "attempts": rows, "errors": errors, "complete_capture": not errors and not pending
              and metadata.get("status") in TERMINAL_WORKFLOW_STATES}
    def encode(value):
        if isinstance(value, Decimal):
            return float(value)
        raise TypeError(type(value).__name__)
    (output / "profiles.json").write_text(json.dumps(result, indent=2, default=encode) + "\n")
    write_table(rows, output / "task_profiles.tsv")
    manifest = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  ./{path.relative_to(output)}"
                for path in sorted(output.rglob("*")) if path.is_file()]
    (output / "evidence-manifest.sha256").write_text("\n".join(manifest) + "\n")
    unexpected = [error for error in errors if not (error.get("expected_after_preemption") or error.get("expected_after_cache_hit")
                                                    or error.get("expected_after_unfinished_attempt"))]
    if unexpected:
        raise ValueError(f"incomplete capture: {len(unexpected)} unexpected unavailable artifacts; see profiles.json")
    return result


def write_table(rows, path):
    """One line per attempt: resources requested, peaks observed and input sizes."""
    peaks = ("mean_cores", "peak_working_set_gib", "peak_memory_anon_gib", "peak_disk_used_gib")
    columns = ["call", "shard", "attempt", "sample", "status", "cache_hit", "preemptible", "machine_type",
               "cpu", "memory", "disks", "start", "end", "reference_release", "input_gib", "input_gib_by_name", *peaks]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            runtime, monitoring = row["runtime"] or {}, row["monitoring"] or {}
            values = [row["call"], row["shard"], row["attempt"], row["sample"], row["status"],
                      bool((row["cache"] or {}).get("hit")), runtime.get("preemptible"),
                      runtime.get("predefinedMachineType"), runtime.get("cpu"), runtime.get("memory"),
                      runtime.get("disks"), row["start"], row["end"], row["reference_release"],
                      round(sum(row["input_bytes"].values()) / 2**30, 4),
                      json.dumps({key: round(value / 2**30, 4) for key, value in row["input_bytes"].items()}, sort_keys=True)] + \
                     [None if monitoring.get(key) is None else round(float(monitoring[key]), 4) for key in peaks]
            handle.write("\t".join("" if value is None else str(value) for value in values) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("output", type=Path, help="new directory for this capture")
    parser.add_argument("--skip-batch", action="store_true", help="omit Batch job descriptions (needs batch.jobs.get)")
    parser.add_argument("--workers", type=int, default=8, help="maximum concurrent cloud reads (default: 8)")
    parser.add_argument("--transport", choices=("rest", "gcloud"), default="rest", help="cloud read transport")
    args = parser.parse_args()
    if args.transport == "rest":
        from profile_cloud import ProfileCloud
        with ProfileCloud() as reader:
            result = collect(args.metadata, args.output, workers=args.workers, batch=not args.skip_batch, cloud_read=reader)
            print(f"Cloud requests: {reader.requests}; reused reads: {reader.cache_hits}")
    else:
        result = collect(args.metadata, args.output, workers=args.workers, batch=not args.skip_batch)
    print(f"Captured {len(result['attempts'])} terminal attempts; {len(result['not_terminal'])} still pending")
    if result["errors"]:
        print(f"WARNING: {len(result['errors'])} logs of cache hits or unfinished attempts unavailable; capture remains incomplete")
