#!/usr/bin/env python3
"""Summarize captured Cromwell/Batch attempts with frozen GCP list prices."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RATES = REPO_ROOT / "config/backends/gcp/gcp-rates-americas-20260830.json"
CUSTOM_N2 = re.compile(r"^n2-custom-(\d+)-(\d+)(?:-ext)?$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
ZERO = Decimal(0)
HOUR = Decimal(3600)
GIB = Decimal(1024) ** 3
KIB_PER_GIB = Decimal(1024) ** 2
PIPELINE_PHASES = {
    "pretrim_fastqc": "fastq_qc",
    "posttrim_fastqc": "fastq_qc",
    "mqc": "fastq_qc",
    "aumi": "read_preparation",
    "cutadapt_umi": "read_preparation",
    "cutadapt_noumi": "read_preparation",
    "star_align": "alignment",
    "feature_counts": "conventional_expression",
    "rsem_quant": "conventional_expression",
    "combined_contamination_qc": "contamination_qc",
    "bowtie2_globin": "contamination_qc",
    "bowtie2_rrna": "contamination_qc",
    "bowtie2_phix": "contamination_qc",
    "md": "alignment_qc",
    "rnaqc": "alignment_qc",
    "chrinfo": "alignment_qc",
    "udup": "umi_molecule_expression",
    "umi_molecule_feature_counts_task": "umi_molecule_expression",
    "umi_molecule_rsem": "umi_molecule_expression",
    "qc_report": "reporting_gather",
    "merge_results": "reporting_gather",
    "merge_umi_expression": "reporting_gather",
    "mqc_pa": "reporting_gather",
}
PIPELINE_PHASE_LABELS = {
    "fastq_qc": "FASTQ QC",
    "read_preparation": "Read preparation",
    "alignment": "Alignment",
    "conventional_expression": "Conventional expression",
    "contamination_qc": "Contamination QC",
    "alignment_qc": "Alignment QC",
    "umi_molecule_expression": "UMI processing / molecule expression",
    "reporting_gather": "Reporting / gather",
}


def fail(message: str) -> None:
    raise ValueError(message)


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        fail(f"expected JSON object: {path}")
    return value


def dec(value: object, label: str, positive: bool = False) -> Decimal:
    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError):
        fail(f"invalid {label}: {value!r}")
    if not number.is_finite() or number < 0 or (positive and number <= 0):
        fail(f"invalid {label}: {value!r}")
    return number


def out(value: Decimal | None) -> float | None:
    return None if value is None else float(value.quantize(Decimal("0.000000000001")))


def timestamp(value: object, label: str) -> datetime:
    if not isinstance(value, str):
        fail(f"missing {label}")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        fail(f"invalid {label}: {value!r}")
    if parsed.tzinfo is None:
        fail(f"timestamp lacks a time zone for {label}")
    return parsed


def elapsed(start: datetime, end: datetime, label: str) -> Decimal:
    seconds = Decimal(str((end - start).total_seconds()))
    if seconds < 0:
        fail(f"negative duration for {label}")
    return seconds


def batch_duration(value: object, label: str) -> Decimal:
    if not isinstance(value, str) or not value.endswith("s"):
        fail(f"invalid {label}: {value!r}")
    return dec(value[:-1], label, positive=True)


def verify_evidence(evidence: Path) -> None:
    manifest = evidence / "evidence-manifest.sha256"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    evidence_root = evidence.resolve()
    listed = set()
    for line in lines:
        parts = line.split("  ", 1)
        if len(parts) != 2 or not SHA256.fullmatch(parts[0]) or not parts[1].startswith("./"):
            fail(f"invalid evidence checksum record: {line!r}")
        relative = Path(parts[1][2:])
        artifact = (evidence / relative).resolve()
        if evidence_root not in artifact.parents:
            fail(f"invalid evidence artifact path: {relative}")
        if not artifact.is_file():
            fail(f"missing evidence artifact: {relative}")
        if relative in listed or hashlib.sha256(artifact.read_bytes()).hexdigest() != parts[0]:
            fail(f"evidence checksum mismatch: {relative}")
        listed.add(relative)
    present = {
        path.relative_to(evidence)
        for path in evidence.rglob("*")
        if path.is_file() and path != manifest
    }
    if not listed or listed != present:
        fail("evidence checksum manifest does not match captured files")


def load_rates(path: Path) -> tuple[dict, dict[str, dict[str, Decimal]], dict[str, Decimal]]:
    rates = load_json(path)
    source = rates.get("source", {})
    if (
        rates.get("schema_version") != 1
        or rates.get("currency") != "USD"
        or rates.get("machine_family") != "N2 custom"
        or not SHA256.fullmatch(str(source.get("snapshot_sha256", "")))
    ):
        fail("unsupported or unprovenanced GCP rate manifest")
    compute = {}
    for market in ("STANDARD", "SPOT"):
        values = rates.get("compute", {}).get(market, {})
        compute[market] = {
            "vcpu": dec(values.get("vcpu_hour_usd"), f"{market} vCPU rate", True),
            "memory": dec(values.get("memory_gib_hour_usd"), f"{market} memory rate", True),
        }
    disk_section = rates.get("disk", {})
    month_hours = dec(disk_section.get("hours_per_month"), "month hours", True)
    disk = {}
    for kind in ("pd-standard", "pd-balanced", "pd-ssd"):
        disk[kind] = dec(
            disk_section.get(kind, {}).get("gib_month_usd"), f"{kind} rate", True
        ) / month_hours
    return rates, compute, disk


def only(value: object, label: str) -> dict:
    if not isinstance(value, list) or len(value) != 1 or not isinstance(value[0], dict):
        fail(f"expected exactly one {label}")
    return value[0]


def phase_times(job: dict, terminal: str) -> tuple[Decimal, Decimal]:
    transitions = {}
    for event in job.get("status", {}).get("statusEvents", []):
        match = re.search(
            r"Job state is set from ([A-Z_]+) to ([A-Z_]+)",
            str(event.get("description", "")),
        )
        if match:
            key = (match.group(1), match.group(2))
            if key in transitions:
                fail(f"{job.get('name')}: duplicate Batch transition {key}")
            transitions[key] = timestamp(event.get("eventTime"), str(key))
    try:
        scheduled = transitions[("QUEUED", "SCHEDULED")]
        running = transitions[("SCHEDULED", "RUNNING")]
        transitions[("RUNNING", terminal)]
    except KeyError:
        fail(f"{job.get('name')}: incomplete Batch status transitions")
    created = timestamp(job.get("createTime"), "Batch createTime")
    return elapsed(created, scheduled, "Batch queue"), elapsed(
        scheduled, running, "Batch provisioning"
    )


def monitor(path: Path) -> dict[str, Decimal | int | None]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "epoch_s",
            "cpu_usage_usec",
            "memory_current_bytes",
            "memory_peak_bytes",
            "disk_used_kb",
            "disk_available_kb",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fail(f"unexpected monitoring header: {path}")
        rows = list(reader)
    if not rows:
        fail(f"monitoring file contains no samples: {path}")

    def values(key: str) -> list[Decimal]:
        return [
            dec(row[key], f"monitoring {key}")
            for row in rows
            if row.get(key) not in {None, "", "NA", "max"}
        ]

    epochs = values("epoch_s")
    if len(epochs) != len(rows) or epochs != sorted(epochs):
        fail(f"missing or unordered monitoring epochs: {path}")
    cpu = values("cpu_usage_usec")
    if cpu and (len(cpu) != len(rows) or cpu != sorted(cpu)):
        fail(f"missing or unordered monitoring CPU counters: {path}")
    observed = epochs[-1] - epochs[0]
    cpu_seconds = (cpu[-1] - cpu[0]) / Decimal(1_000_000) if cpu else None
    memory = values("memory_current_bytes") + values("memory_peak_bytes")
    disk_used = values("disk_used_kb")
    disk_free = values("disk_available_kb")
    return {
        "sample_count": len(rows),
        "observed_seconds": observed,
        "cpu_seconds": cpu_seconds,
        "mean_cores": cpu_seconds / observed if cpu_seconds is not None and observed else None,
        "peak_memory_gib": max(memory) / GIB if memory else None,
        "peak_disk_used_gib": max(disk_used) / KIB_PER_GIB if disk_used else None,
        "minimum_disk_free_gib": min(disk_free) / KIB_PER_GIB if disk_free else None,
    }


def summarize_attempt(
    evidence: Path,
    call: str,
    attempt: dict,
    jobs: dict[str, dict],
    compute_rates: dict[str, dict[str, Decimal]],
    disk_rates: dict[str, Decimal],
) -> dict:
    attempt_number = attempt.get("attempt")
    shard = attempt.get("shardIndex")
    if (
        not re.fullmatch(r"[A-Za-z0-9_.-]+", call)
        or not isinstance(attempt_number, int)
        or attempt_number < 1
        or not isinstance(shard, int)
    ):
        fail(f"invalid call attempt: {call}")
    job_id = attempt.get("jobId")
    if not isinstance(job_id, str) or job_id not in jobs:
        fail(f"{call} attempt {attempt_number}: missing captured Batch job")
    job = jobs[job_id]
    status = job.get("status", {})
    terminal = status.get("state")
    if terminal not in {"SUCCEEDED", "FAILED"}:
        fail(f"{job_id}: Batch job is not terminal")
    queue_seconds, provisioning_seconds = phase_times(job, terminal)
    running_seconds = batch_duration(status.get("runDuration"), "Batch runDuration")

    actual = only(
        status.get("taskGroups", {}).get("group0", {}).get("instances"),
        "actual Batch instance",
    )
    policy = only(
        job.get("allocationPolicy", {}).get("instances"), "Batch instance policy"
    ).get("policy", {})
    machine = actual.get("machineType")
    match = CUSTOM_N2.fullmatch(str(machine))
    if not match or machine != policy.get("machineType"):
        fail(f"{job_id}: unsupported or inconsistent machine type")
    vcpu = int(match.group(1))
    memory_gib = Decimal(match.group(2)) / Decimal(1024)
    market = actual.get("provisioningModel")
    if market not in compute_rates or market != policy.get("provisioningModel"):
        fail(f"{job_id}: unsupported or inconsistent provisioning model")

    boot = actual.get("bootDisk", {})
    disks = [
        {
            "role": "boot",
            "type": boot.get("type"),
            "size_gib": int(dec(boot.get("sizeGb"), "boot disk size", True)),
        }
    ]
    for allocation in policy.get("disks", []):
        new_disk = allocation.get("newDisk", {})
        disks.append(
            {
                "role": allocation.get("deviceName", "work"),
                "type": new_disk.get("type"),
                "size_gib": int(dec(new_disk.get("sizeGb"), "work disk size", True)),
            }
        )
    if any(item["type"] not in disk_rates for item in disks):
        fail(f"{job_id}: unpriced disk type")

    short_call = call.removeprefix("rnaseq_pipeline.")
    pipeline_phase = PIPELINE_PHASES.get(short_call)
    if pipeline_phase is None:
        fail(f"{call}: no fixed pipeline-phase assignment")
    monitor_path = evidence / "task-streams" / (
        f"{short_call}.shard-{shard}.attempt-{attempt_number}.monitoring"
    )
    metrics = monitor(monitor_path)
    observed = min(metrics["observed_seconds"], running_seconds)
    outside_monitor = running_seconds - observed
    hours = running_seconds / HOUR
    vcpu_cost = Decimal(vcpu) * hours * compute_rates[market]["vcpu"]
    memory_cost = memory_gib * hours * compute_rates[market]["memory"]
    disk_gib = sum(Decimal(item["size_gib"]) for item in disks)
    disk_cost = sum(
        Decimal(item["size_gib"]) * hours * disk_rates[item["type"]]
        for item in disks
    )
    cost = vcpu_cost + memory_cost + disk_cost
    failed = terminal != "SUCCEEDED" or attempt.get("executionStatus") != "Done"
    metadata_seconds = elapsed(
        timestamp(attempt.get("start"), f"{call} start"),
        timestamp(attempt.get("end"), f"{call} end"),
        f"{call} metadata",
    )
    return {
        "call": call,
        "pipeline_phase": pipeline_phase,
        "shard": shard,
        "attempt": attempt_number,
        "execution_status": attempt.get("executionStatus"),
        "batch_state": terminal,
        "failed_work": failed,
        "batch_job": job_id,
        "provisioning_model": market,
        "machine_type": machine,
        "provisioned_vcpu": vcpu,
        "provisioned_memory_gib": out(memory_gib),
        "disks": disks,
        "phase_seconds": {
            "metadata_elapsed": out(metadata_seconds),
            "batch_queue": out(queue_seconds),
            "batch_provisioning": out(provisioning_seconds),
            "batch_running": out(running_seconds),
            "monitor_observed": out(observed),
            "running_outside_monitor_window": out(outside_monitor),
        },
        "monitoring": {
            key: out(value) if isinstance(value, Decimal) else value
            for key, value in metrics.items()
        },
        "modeled_worker_cost_usd": {
            "vcpu": out(vcpu_cost),
            "memory": out(memory_cost),
            "disk": out(disk_cost),
            "total": out(cost),
        },
        "_raw": {
            "queue": queue_seconds,
            "provisioning": provisioning_seconds,
            "running": running_seconds,
            "observed": observed,
            "outside": outside_monitor,
            "vcpu_hours": Decimal(vcpu) * hours,
            "memory_hours": memory_gib * hours,
            "disk_hours": disk_gib * hours,
            "vcpu_cost": vcpu_cost,
            "memory_cost": memory_cost,
            "disk_cost": disk_cost,
            "cost": cost,
            "failed_cost": cost if failed else ZERO,
        },
    }


def rollup(attempts: list[dict]) -> dict:
    totals = defaultdict(lambda: ZERO)
    models = defaultdict(int)
    for attempt in attempts:
        for key, value in attempt["_raw"].items():
            totals[key] += value
        models[attempt["provisioning_model"]] += 1
    return {
        "attempt_count": len(attempts),
        "failed_attempt_count": sum(item["failed_work"] for item in attempts),
        "provisioning_model_counts": dict(sorted(models.items())),
        "batch_running_seconds": out(totals["running"]),
        "provisioned_vcpu_hours": out(totals["vcpu_hours"]),
        "provisioned_memory_gib_hours": out(totals["memory_hours"]),
        "provisioned_disk_gib_hours": out(totals["disk_hours"]),
        "modeled_worker_cost_usd": {
            "vcpu": out(totals["vcpu_cost"]),
            "memory": out(totals["memory_cost"]),
            "disk": out(totals["disk_cost"]),
            "total": out(totals["cost"]),
            "failed_work": out(totals["failed_cost"]),
        },
    }


def summarize(evidence: Path, rates_path: Path) -> dict:
    verify_evidence(evidence)
    metadata = load_json(evidence / "metadata.json")
    capture = load_json(evidence / "capture-status.json")
    repository = load_json(evidence / "repository.json")
    workflow_id = metadata.get("id")
    if (
        capture.get("complete") is not True
        or capture.get("workflow_id") != workflow_id
        or capture.get("repository_revision") != repository.get("revision")
        or repository.get("expected_submission_revision") != repository.get("revision")
        or capture.get("repository_clean") != repository.get("clean")
    ):
        fail("capture status is incomplete or inconsistent")

    with (evidence / "input-objects.json").open(encoding="utf-8") as handle:
        input_objects = json.load(handle)
    if not isinstance(input_objects, list) or len(input_objects) != capture.get(
        "submitted_gcs_object_count"
    ):
        fail("input-object manifest is inconsistent with capture status")
    submitted_bytes = sum(
        dec(item.get("size_bytes"), "input object size") for item in input_objects
    )
    with (evidence / "output-objects.json").open(encoding="utf-8") as handle:
        output_objects = json.load(handle)
    if not isinstance(output_objects, list) or len(output_objects) != capture.get(
        "top_level_output_object_count"
    ):
        fail("output-object manifest is inconsistent with capture status")
    output_bytes = sum(
        dec(item.get("size_bytes"), "output object size") for item in output_objects
    )

    rates, compute_rates, disk_rates = load_rates(rates_path)
    jobs = {}
    for path in sorted((evidence / "batch-jobs").glob("*.json")):
        job = load_json(path)
        name = job.get("name")
        if not isinstance(name, str) or name in jobs:
            fail(f"invalid or duplicate Batch job name: {path}")
        jobs[name] = job

    attempts = []
    calls = metadata.get("calls")
    if not isinstance(calls, dict):
        fail("metadata lacks calls")
    for call in sorted(calls):
        if not isinstance(calls[call], list):
            fail(f"metadata attempts are invalid for {call}")
        for attempt in sorted(
            calls[call], key=lambda item: (item.get("shardIndex"), item.get("attempt"))
        ):
            attempts.append(
                summarize_attempt(
                    evidence, call, attempt, jobs, compute_rates, disk_rates
                )
            )
    if (
        not attempts
        or len(attempts) != capture.get("attempt_count")
        or len(jobs) != len(attempts)
    ):
        fail("captured attempts and Batch jobs are inconsistent")

    totals = defaultdict(lambda: ZERO)
    for attempt in attempts:
        for key, value in attempt["_raw"].items():
            totals[key] += value
    by_call = {
        call: rollup([item for item in attempts if item["call"] == call])
        for call in sorted(calls)
    }
    by_market = {
        market: rollup(
            [item for item in attempts if item["provisioning_model"] == market]
        )
        for market in sorted({item["provisioning_model"] for item in attempts})
    }
    by_pipeline_phase = {
        phase: {
            "label": PIPELINE_PHASE_LABELS[phase],
            **rollup([item for item in attempts if item["pipeline_phase"] == phase]),
        }
        for phase in PIPELINE_PHASE_LABELS
        if any(item["pipeline_phase"] == phase for item in attempts)
    }
    total_rollup = rollup(attempts)
    observed_cost = sum(
        item["_raw"]["cost"] * item["_raw"]["observed"] / item["_raw"]["running"]
        for item in attempts
    )
    for attempt in attempts:
        del attempt["_raw"]
    total_cost = totals["cost"]
    return {
        "schema_version": 1,
        "workflow": {
            "id": workflow_id,
            "status": metadata.get("status"),
            "wall_seconds": out(
                elapsed(
                    timestamp(metadata.get("start"), "workflow start"),
                    timestamp(metadata.get("end"), "workflow end"),
                    "workflow",
                )
            ),
            "repository_revision": repository.get("revision"),
            "repository_clean": repository.get("clean"),
            "submitted_gcs_object_count": len(input_objects),
            "submitted_gcs_object_bytes": int(submitted_bytes),
            "top_level_output_object_count": len(output_objects),
            "top_level_output_object_bytes": int(output_bytes),
        },
        "pricing": {
            "currency": rates["currency"],
            "region_scope": rates["region_scope"],
            "machine_family": rates["machine_family"],
            "rate_manifest_sha256": hashlib.sha256(rates_path.read_bytes()).hexdigest(),
            "source": rates["source"],
        },
        "totals": total_rollup,
        "by_call": by_call,
        "by_market": by_market,
        "by_pipeline_phase": by_pipeline_phase,
        "phase_totals": {
            "batch_queue": {"seconds": out(totals["queue"]), "worker_cost_usd": None},
            "batch_provisioning": {
                "seconds": out(totals["provisioning"]),
                "worker_cost_usd": None,
            },
            "batch_running": {
                "seconds": out(totals["running"]),
                "worker_cost_usd": out(total_cost),
                "monitor_observed": {
                    "seconds": out(totals["observed"]),
                    "cost_equivalent_usd": out(observed_cost),
                },
                "outside_monitor_window": {
                    "seconds": out(totals["outside"]),
                    "cost_equivalent_usd": out(total_cost - observed_cost),
                },
            },
        },
        "attempts": attempts,
        "cost_scope": {
            "basis": (
                "Actual Batch runDuration, actual N2 custom shape/disks, and "
                "frozen public list rates."
            ),
            "failed_spot_work_included": True,
            "excluded": [
                "controller VM, Cloud Storage, logging/monitoring, and network",
                "taxes, credits, discounts, and the account-level disk free tier",
                "minimum-duration billing adjustments",
            ],
            "phase_caveat": (
                "The monitoring window is only temporal coverage; it does not isolate "
                "localization, command, and delocalization."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence_directory", type=Path)
    parser.add_argument("--rates", type=Path, default=DEFAULT_RATES)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = summarize(args.evidence_directory, args.rates)
        rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output:
            with args.output.open("x", encoding="utf-8") as handle:
                handle.write(rendered)
        else:
            sys.stdout.write(rendered)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
