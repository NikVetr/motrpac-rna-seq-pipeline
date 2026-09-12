"""Require regional execution storage; report remote source buckets before launch."""
import argparse
import json
import os
import subprocess
from pathlib import Path
from urllib.parse import urlsplit


def buckets(value):
    if isinstance(value, dict):
        return set().union(*(buckets(item) for item in value.values()))
    if isinstance(value, list):
        return set().union(*(buckets(item) for item in value))
    if isinstance(value, str) and value.startswith("gs://"):
        bucket = urlsplit(value).netloc
        if not bucket:
            raise ValueError("GCS URI lacks a bucket")
        return {bucket}
    return set()


def check(root, region, options, inputs, describe):
    roots = buckets(root)
    if len(roots) != 1 or not region:
        raise ValueError("set an explicit GCS execution root and Batch region")
    override = options.get("gcp_batch_gcs_root", root)
    if override.rstrip("/") != root.rstrip("/"):
        raise ValueError("workflow gcp_batch_gcs_root differs from the checked execution root")
    zones = options.get("default_runtime_attributes", {}).get("zones", "")
    if not isinstance(zones, str) or not zones.split():
        raise ValueError("workflow options must specify worker zones explicitly")
    if any(zone.rsplit("-", 1)[0] != region for zone in zones.split()):
        raise ValueError("worker zones differ from the Batch region")
    if inputs.get("rnaseq_pipeline.use_e2", False):
        runtime = options.get("default_runtime_attributes", {})
        if inputs.get("rnaseq_pipeline.prefer_predefined_n1", False) or runtime.get("cpuPlatform") or runtime.get("predefinedMachineType"):
            raise ValueError("E2 policy requires predefined N1 disabled and no CPU platform/machine override")
    if inputs.get("rnaseq_pipeline.prefer_predefined_n1", False):
        runtime = options.get("default_runtime_attributes", {})
        if region != "us-west2" or runtime.get("cpuPlatform") or runtime.get("predefinedMachineType"):
            raise ValueError("predefined N1 policy requires us-west2 and no CPU platform/machine override")
    report = []
    for bucket in sorted(roots | buckets(inputs) | buckets(options)):
        obj = describe(bucket)
        location = obj["location"].lower()
        local = location == region and obj["locationType"] == "region"
        execution = bucket in roots
        if execution and not local:
            raise ValueError(f"execution bucket {bucket} is {location}; require regional {region}")
        report.append({"bucket": bucket, "location": location,
                       "role": "execution" if execution else "source", "colocated": local})
    return report


def describe_bucket(bucket):
    return json.loads(subprocess.check_output(
        ["gcloud", "storage", "buckets", "describe", "gs://" + bucket,
         "--raw", "--format=json"], text=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--options", type=Path, required=True)
    parser.add_argument("--execution-root", default=os.environ.get("MOTRPAC_GCP_BATCH_ROOT"))
    parser.add_argument("--region", default=os.environ.get("MOTRPAC_GCP_BATCH_LOCATION"))
    parser.add_argument("command", nargs=argparse.REMAINDER, help="optional launch command after --")
    args = parser.parse_args()
    report = check(args.execution_root, args.region, json.loads(args.options.read_text()),
                   json.loads(args.inputs.read_text()), describe_bucket)
    print(json.dumps(report, indent=2), flush=True)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if command:
        # The checked values must also reach the Cromwell backend.
        env = dict(os.environ, MOTRPAC_GCP_BATCH_ROOT=args.execution_root,
                   MOTRPAC_GCP_BATCH_LOCATION=args.region)
        return subprocess.call(command, env=env)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
