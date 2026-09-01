import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/gcp/summarize_workflow_cost.py"
RATES = REPO_ROOT / "config/backends/gcp/gcp-rates-americas-20260830.json"


class GcpCostSummaryTests(unittest.TestCase):
    def make_evidence(self, root: Path) -> Path:
        evidence = root / "evidence"
        (evidence / "batch-jobs").mkdir(parents=True)
        (evidence / "task-streams").mkdir()
        workflow_id = "12345678-1234-1234-1234-123456789abc"
        revision = "a" * 40
        attempts = []
        cases = (
            (1, "Failed", "FAILED", "SPOT", 100),
            (2, "Done", "SUCCEEDED", "STANDARD", 200),
        )
        for attempt, execution, state, market, seconds in cases:
            job_name = f"job-test-{attempt}"
            job_id = f"projects/test/locations/us-central1/jobs/{job_name}"
            start_second = attempt * 10
            attempts.append(
                {
                    "attempt": attempt,
                    "executionStatus": execution,
                    "jobId": job_id,
                    "shardIndex": 0,
                    "start": f"2026-08-30T00:00:{start_second:02d}Z",
                    "end": "2026-08-30T00:05:00Z",
                }
            )
            terminal_time = "2026-08-30T00:01:45Z" if seconds == 100 else "2026-08-30T00:03:25Z"
            job = {
                "name": job_id,
                "createTime": "2026-08-30T00:00:00Z",
                "status": {
                    "runDuration": f"{seconds}s",
                    "state": state,
                    "statusEvents": [
                        {
                            "description": "Job state is set from QUEUED to SCHEDULED",
                            "eventTime": "2026-08-30T00:00:02Z",
                        },
                        {
                            "description": "Job state is set from SCHEDULED to RUNNING",
                            "eventTime": "2026-08-30T00:00:05Z",
                        },
                        {
                            "description": f"Job state is set from RUNNING to {state}",
                            "eventTime": terminal_time,
                        },
                    ],
                    "taskGroups": {
                        "group0": {
                            "instances": [
                                {
                                    "bootDisk": {
                                        "sizeGb": "41",
                                        "type": "pd-balanced",
                                    },
                                    "machineType": "n2-custom-2-8192",
                                    "provisioningModel": market,
                                }
                            ]
                        }
                    },
                },
                "allocationPolicy": {
                    "instances": [
                        {
                            "policy": {
                                "disks": [
                                    {
                                        "deviceName": "local-disk",
                                        "newDisk": {
                                            "sizeGb": "120",
                                            "type": "pd-ssd",
                                        },
                                    }
                                ],
                                "machineType": "n2-custom-2-8192",
                                "provisioningModel": market,
                            }
                        }
                    ]
                },
            }
            (evidence / "batch-jobs" / f"{job_name}.json").write_text(
                json.dumps(job), encoding="utf-8"
            )
            monitor = evidence / "task-streams" / (
                f"star_align.shard-0.attempt-{attempt}.monitoring"
            )
            monitor.write_text(
                "timestamp_utc\tepoch_s\tcpu_usage_usec\tmemory_current_bytes\t"
                "memory_peak_bytes\tmemory_limit_bytes\thost_mem_available_kb\t"
                "disk_used_kb\tdisk_available_kb\n"
                "2026-08-30T00:00:10Z\t10\t1000000\t1073741824\t1073741824\t"
                "max\t1\t1048576\t2097152\n"
                "2026-08-30T00:00:40Z\t40\t31000000\t2147483648\t2147483648\t"
                "max\t1\t2097152\t1048576\n",
                encoding="utf-8",
            )

        metadata = {
            "id": workflow_id,
            "status": "Succeeded",
            "start": "2026-08-30T00:00:00Z",
            "end": "2026-08-30T00:06:00Z",
            "calls": {"rnaseq_pipeline.star_align": attempts},
        }
        (evidence / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        (evidence / "repository.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "revision": revision,
                    "expected_submission_revision": revision,
                    "clean": False,
                }
            ),
            encoding="utf-8",
        )
        (evidence / "capture-status.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "workflow_id": workflow_id,
                    "repository_revision": revision,
                    "repository_clean": False,
                    "attempt_count": 2,
                    "submitted_gcs_object_count": 2,
                    "top_level_output_object_count": 2,
                    "missing_artifact_count": 0,
                    "complete": True,
                }
            ),
            encoding="utf-8",
        )
        (evidence / "input-objects.json").write_text(
            json.dumps(
                [
                    {"uri": "gs://test/R1.fastq.gz", "size_bytes": "1000"},
                    {"uri": "gs://test/R2.fastq.gz", "size_bytes": "2000"},
                ]
            ),
            encoding="utf-8",
        )
        (evidence / "output-objects.json").write_text(
            json.dumps(
                [
                    {"uri": "gs://test/counts.txt", "size_bytes": "3000"},
                    {"uri": "gs://test/qc.csv", "size_bytes": "4000"},
                ]
            ),
            encoding="utf-8",
        )
        manifest_lines = []
        for path in sorted(item for item in evidence.rglob("*") if item.is_file()):
            manifest_lines.append(
                "{}  ./{}".format(
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    path.relative_to(evidence),
                )
            )
        (evidence / "evidence-manifest.sha256").write_text(
            "\n".join(manifest_lines) + "\n", encoding="utf-8"
        )
        return evidence

    def test_includes_failed_spot_attempt_and_phase_costs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            evidence = self.make_evidence(root)
            output = root / "summary.json"
            subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence), "--output", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))

            self.assertEqual(2, summary["totals"]["attempt_count"])
            self.assertEqual(1, summary["totals"]["failed_attempt_count"])
            self.assertEqual(
                2, summary["by_pipeline_phase"]["alignment"]["attempt_count"]
            )
            self.assertEqual(300, summary["phase_totals"]["batch_running"]["seconds"])
            self.assertEqual(
                60,
                summary["phase_totals"]["batch_running"]["monitor_observed"]["seconds"],
            )
            self.assertEqual(3000, summary["workflow"]["submitted_gcs_object_bytes"])
            self.assertEqual(7000, summary["workflow"]["top_level_output_object_bytes"])
            self.assertTrue(summary["cost_scope"]["failed_spot_work_included"])
            self.assertTrue(summary["attempts"][0]["failed_work"])
            self.assertEqual("SPOT", summary["attempts"][0]["provisioning_model"])
            self.assertEqual(
                summary["by_market"]["SPOT"]["modeled_worker_cost_usd"]["total"],
                summary["totals"]["modeled_worker_cost_usd"]["failed_work"],
            )
            self.assertGreater(
                summary["totals"]["modeled_worker_cost_usd"]["total"],
                summary["totals"]["modeled_worker_cost_usd"]["failed_work"],
            )
            self.assertEqual(
                hashlib.sha256(RATES.read_bytes()).hexdigest(),
                summary["pricing"]["rate_manifest_sha256"],
            )

    def test_tampered_evidence_fails_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            evidence = self.make_evidence(Path(temporary_directory))
            (evidence / "batch-jobs" / "job-test-1.json").unlink()
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(0, result.returncode)
            self.assertIn(
                "missing evidence artifact",
                result.stderr,
            )


if __name__ == "__main__":
    unittest.main()
