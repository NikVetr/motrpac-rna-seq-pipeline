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
from functools import partial
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import capture_task_profiles  # noqa: E402
from make_merge_inputs import workflow_metadata  # noqa: E402
from summarize_workflow_cost import verify_evidence  # noqa: E402

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


def reuse_capture(root, metadata, target, sources, batch):
    """Reuse explicitly supplied terminal captures only after checking metadata and all hashes."""
    workflow = root.rsplit("/", 1)[-1]
    for source in sources:
        directory = source / workflow
        if not directory.is_dir():
            continue
        verify_evidence(directory)
        previous = json.loads((directory / "metadata.json").read_text())
        if previous != metadata:
            print(f"Not reusing {directory}: workflow metadata differs; capturing again", file=sys.stderr)
            continue
        result = json.loads((directory / "profiles.json").read_text())
        if (metadata.get("id") != workflow or metadata.get("workflowRoot", "").rstrip("/") != root or
                result.get("workflow_id") != workflow or result.get("workflow_status") != metadata.get("status") or
                metadata.get("status") not in capture_task_profiles.TERMINAL_WORKFLOW_STATES or result.get("not_terminal")):
            raise ValueError(f"invalid or nonterminal reusable capture: {directory}")
        if any(not (error.get("expected_after_preemption") or error.get("expected_after_cache_hit") or
                    error.get("expected_after_unfinished_attempt")) for error in result["errors"]):
            raise ValueError(f"reusable capture has unexpected missing evidence: {directory}")
        attempts = {(row["call"], row["shard"], row["attempt"]): row for row in result["attempts"]}
        expected = {(call, a["shardIndex"], a["attempt"]): a
                    for call, calls in metadata["calls"].items() for a in calls}
        if not expected or len(attempts) != len(result["attempts"]) or set(attempts) != set(expected):
            raise ValueError(f"reusable capture lacks attempts: {directory}")
        if any(a.get("executionStatus") not in capture_task_profiles.TERMINAL_ATTEMPT_STATES or
               attempts[key].get("status") != a["executionStatus"] for key, a in expected.items()):
            raise ValueError(f"reusable capture has inconsistent attempt states: {directory}")
        if batch and any(a.get("jobId") and not attempts[key].get("batch") for key, a in expected.items()):
            print(f"Not reusing {directory}: Batch descriptions missing; capturing again", file=sys.stderr)
            continue
        shutil.copytree(directory, target)
        verify_evidence(target)
        print(f"Reused verified capture: {directory}")
        return result, str(directory)
    return None


def capture(roots, output, batch=True, read_metadata=workflow_metadata, collect=capture_task_profiles.collect, reuse=()):
    roots = list(dict.fromkeys(root.rstrip("/") for root in roots))
    if not roots:
        raise ValueError("no workflow run roots to capture")
    output.mkdir(parents=True, exist_ok=False)
    failures, reused = {}, {}
    for root in roots:  # grows as cache origins are found
        try:
            metadata = read_metadata(root)
            if metadata.get("status") not in capture_task_profiles.TERMINAL_WORKFLOW_STATES:
                raise ValueError(f"workflow is not terminal: {root} ({metadata.get('status')})")
            roots.extend(origin for origin in cache_origins(metadata, root) if origin not in roots)
            target = output / root.rsplit("/", 1)[-1]
            previous = reuse_capture(root, metadata, target, reuse, batch)
            if previous:
                result, reused[root] = previous
            else:
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "metadata.json"
                    path.write_text(json.dumps(metadata))
                    result = collect(path, target, batch=batch)
            if result.get("not_terminal"):
                raise ValueError(f"workflow still has pending attempts: {root}")
            if result["errors"]:
                print(f"WARNING: {root}: {len(result['errors'])} expected missing logs; "
                      "capture remains incomplete (see its profiles.json)", file=sys.stderr)
        except (RuntimeError, ValueError, OSError) as error:
            failures[root] = str(error)
    (output / "runs.json").write_text(json.dumps({"run_roots": roots, "failures": failures, "reused": reused}, indent=2) + "\n")
    if failures:
        raise RuntimeError(f"{len(failures)} of {len(roots)} captures failed; see {output / 'runs.json'}")
    bundle = output.with_name(output.name + ".tar.gz")
    with tarfile.open(bundle, "x:gz", compresslevel=1) as archive:
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
    parser.add_argument("--workers", type=int, default=8, help="maximum concurrent cloud reads (default: 8)")
    parser.add_argument("--transport", choices=("rest", "gcloud"), default="rest", help="cloud read transport")
    parser.add_argument("--reuse-captures", type=Path, action="append", default=[],
                        help="completed local cohort capture directory containing workflow-ID subdirectories; repeatable")
    args = parser.parse_args(argv)
    roots = args.run_root + [root for path in args.merge_inputs for root in merge_run_roots(path)]
    if args.workers < 1:
        parser.error("--workers must be positive")
    if any(not source.is_dir() for source in args.reuse_captures):
        parser.error("--reuse-captures must name existing local capture directories")
    collect = partial(capture_task_profiles.collect, workers=args.workers)
    if args.transport == "rest":
        from profile_cloud import ProfileCloud
        with ProfileCloud() as reader:
            bundle = capture(roots, args.output, batch=not args.skip_batch,
                             collect=partial(collect, cloud_read=reader), reuse=args.reuse_captures)
            print(f"Cloud requests: {reader.requests}; reused reads: {reader.cache_hits}")
    else:
        bundle = capture(roots, args.output, batch=not args.skip_batch, collect=collect, reuse=args.reuse_captures)
    print(f"Wrote {bundle}")
    if args.destination:
        target = args.destination.rstrip("/") + "/" + bundle.name
        subprocess.run(["gcloud", "storage", "cp", "--if-generation-match=0", str(bundle), target], check=True)
        print(f"Copied to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
