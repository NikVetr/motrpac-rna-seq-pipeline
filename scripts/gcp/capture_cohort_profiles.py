#!/usr/bin/env python3
"""Capture task profiles for every workflow behind a cohort and bundle them for storage.

Run on the controller before execution directories are deleted: metadata comes from each run's
metadata.json, or from Caper when that export is missing, and Batch job descriptions need
batch.jobs.get. Workflows are named with --run-root and/or taken from a merge inputs JSON (the
run directories its source files resolve to). Every cache hit's source workflow is added in turn,
since only that workflow measured the task. Each workflow gets its own capture_task_profiles.py
snapshot; the bundle is written only if every capture succeeds.
"""

import argparse
import json
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import capture_task_profiles  # noqa: E402
from make_merge_inputs import workflow_metadata  # noqa: E402

RUN_ROOT = re.compile(r"(gs://.+?/rnaseq_pipeline/[^/]+)/call-")
CACHE_HIT = re.compile(r"Cache Hit: ([0-9a-f-]{36}):")


def merge_run_roots(merge_inputs):
    rows = json.loads(merge_inputs.read_text())["rnaseq_merge.source_rows"][1:]
    return [match[1] for match in (RUN_ROOT.match(uri) for _, _, uri in rows) if match]


def cache_origins(metadata, root):
    """Run roots of the workflows whose outputs this workflow's cache hits reused (same Caper output prefix)."""
    hits = (CACHE_HIT.match((call.get("callCaching") or {}).get("result", ""))
            for calls in metadata.get("calls", {}).values() for call in calls)
    return [root.rsplit("/", 1)[0] + "/" + workflow for workflow in sorted({hit[1] for hit in hits if hit})]


def capture(roots, output, batch=True, read_metadata=workflow_metadata, collect=capture_task_profiles.collect):
    roots = list(dict.fromkeys(root.rstrip("/") for root in roots))
    if not roots:
        raise ValueError("no workflow run roots to capture")
    output.mkdir(parents=True, exist_ok=False)
    failures = {}
    for root in roots:  # grows as cache origins are found
        try:
            metadata = read_metadata(root)
            roots.extend(origin for origin in cache_origins(metadata, root) if origin not in roots)
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "metadata.json"
                path.write_text(json.dumps(metadata))
                collect(path, output / root.rsplit("/", 1)[-1], batch=batch)
        except (RuntimeError, ValueError) as error:
            failures[root] = str(error)
    (output / "runs.json").write_text(json.dumps({"run_roots": roots, "failures": failures}, indent=2) + "\n")
    if failures:
        raise RuntimeError(f"{len(failures)} of {len(roots)} captures failed; see {output / 'runs.json'}")
    bundle = output.with_name(output.name + ".tar.gz")
    with tarfile.open(bundle, "x:gz") as archive:
        archive.add(output, arcname=output.name)
    return bundle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new local directory; the bundle is OUTPUT.tar.gz")
    parser.add_argument("--run-root", action="append", default=[], help="gs://.../rnaseq_pipeline/WORKFLOW_ID")
    parser.add_argument("--merge-inputs", type=Path, action="append", default=[],
                        help="make_merge_inputs.py output whose source runs are captured")
    parser.add_argument("--destination", help="gs:// prefix to copy the bundle to; an existing object is not replaced")
    parser.add_argument("--skip-batch", action="store_true", help="omit Batch job descriptions (needs batch.jobs.get)")
    args = parser.parse_args(argv)
    roots = args.run_root + [root for path in args.merge_inputs for root in merge_run_roots(path)]
    bundle = capture(roots, args.output, batch=not args.skip_batch)
    print(f"Wrote {bundle}")
    if args.destination:
        target = args.destination.rstrip("/") + "/" + bundle.name
        subprocess.run(["gcloud", "storage", "cp", "--if-generation-match=0", str(bundle), target], check=True)
        print(f"Copied to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
