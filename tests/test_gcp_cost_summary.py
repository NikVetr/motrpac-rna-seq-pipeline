import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/gcp/summarize_workflow_cost.py"
N1_RATES = REPO_ROOT / "config/backends/gcp/gcp-rates-n1-americas-20260830.json"
N2_RATES = REPO_ROOT / "config/backends/gcp/gcp-rates-americas-20260830.json"


class GcpCostSummaryTests(unittest.TestCase):
    @staticmethod
    def write_json(path: Path, value) -> None:
        path.write_text(json.dumps(value), encoding="utf-8")

    @staticmethod
    def write_manifest(evidence: Path) -> None:
        manifest = evidence / "evidence-manifest.sha256"
        manifest.unlink(missing_ok=True)
        lines = []
        for path in sorted(item for item in evidence.rglob("*") if item.is_file()):
            lines.append(
                "{}  ./{}".format(
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    path.relative_to(evidence),
                )
            )
        manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def make_evidence(
        self, root: Path, machine_type: str = "custom-2-8192"
    ) -> Path:
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
            events = (
                ("QUEUED to SCHEDULED", "2026-08-30T00:00:02Z"),
                ("SCHEDULED to RUNNING", "2026-08-30T00:00:05Z"),
                (f"RUNNING to {state}", terminal_time),
            )
            status_events = [
                {
                    "description": f"Job state is set from {transition}",
                    "eventTime": time,
                }
                for transition, time in events
            ]
            if attempt == 1:
                status_events[-1]["description"] += (
                    ". Task state is updated from RUNNING to FAILED on "
                    "zones/us-west1-a/instances/1 due to Spot VM preemption "
                    "with exit code 50001."
                )
            allocation = {
                "disks": [{
                    "deviceName": "local-disk",
                    "newDisk": {"sizeGb": "120", "type": "pd-ssd"},
                }],
                "machineType": machine_type,
                "provisioningModel": market,
            }
            status_instance = {
                "bootDisk": {"sizeGb": "41", "type": "pd-balanced"},
                "machineType": allocation["machineType"],
                "provisioningModel": market,
            }
            job = {
                "name": job_id,
                "createTime": "2026-08-30T00:00:00.123456789Z",
                "status": {
                    "runDuration": f"{seconds}s",
                    "state": state,
                    "statusEvents": status_events,
                    "taskGroups": {"group0": {"instances": [status_instance]}},
                },
                "allocationPolicy": {"instances": [{"policy": allocation}]},
            }
            self.write_json(evidence / "batch-jobs" / f"{job_name}.json", job)
            monitor = evidence / "task-streams" / (
                f"star_align.shard-0.attempt-{attempt}.monitoring"
            )
            if attempt == 2:
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
            for stream in ("stdout", "stderr"):
                (evidence / "task-streams" / (
                    f"star_align.shard-0.attempt-{attempt}.{stream}"
                )).write_text("captured\n", encoding="utf-8")

        metadata = {
            "id": workflow_id,
            "status": "Succeeded",
            "start": "2026-08-30T00:00:00Z",
            "end": "2026-08-30T00:06:00Z",
            "calls": {"rnaseq_pipeline.star_align": attempts},
        }
        self.write_json(evidence / "metadata.json", metadata)
        self.write_json(
            evidence / "repository.json",
            {
                "schema_version": 1,
                "revision": revision,
                "expected_submission_revision": revision,
                "clean": False,
            },
        )
        self.write_json(
            evidence / "capture-status.json",
            {
                "schema_version": 1,
                "workflow_id": workflow_id,
                "repository_revision": revision,
                "repository_clean": False,
                "attempt_count": 2,
                "submitted_gcs_object_count": 2,
                "top_level_output_object_count": 2,
                "missing_artifact_count": 0,
                "expected_unavailable_artifact_count": 1,
                "complete": True,
            },
        )
        self.write_json(
            evidence / "input-objects.json",
            [
                {"uri": "gs://test/R1.fastq.gz", "size_bytes": "1000"},
                {"uri": "gs://test/R2.fastq.gz", "size_bytes": "2000"},
            ],
        )
        self.write_json(
            evidence / "output-objects.json",
            [
                {"uri": "gs://test/counts.txt", "size_bytes": "3000"},
                {"uri": "gs://test/qc.csv", "size_bytes": "4000"},
            ],
        )
        (evidence / "expected-unavailable-artifacts.tsv").write_text(
            "kind\tcall\tshard\tattempt\turi\treason\n"
            "monitoring\trnaseq_pipeline.star_align\t0\t1\t"
            "gs://test/job-test-1/monitoring.log\t"
            "recognized_batch_infrastructure_failure\n",
            encoding="utf-8",
        )
        self.write_manifest(evidence)
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
                30,
                summary["phase_totals"]["batch_running"]["monitor_observed"]["seconds"],
            )
            self.assertEqual(3000, summary["workflow"]["submitted_gcs_object_bytes"])
            self.assertEqual(7000, summary["workflow"]["top_level_output_object_bytes"])
            self.assertTrue(summary["cost_scope"]["failed_spot_work_included"])
            self.assertTrue(summary["attempts"][0]["failed_work"])
            self.assertIsNone(summary["attempts"][0]["monitoring"])
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
                "N1 custom",
                summary["pricing"]["machine_family"],
            )
            self.assertIn("actual N1 custom", summary["cost_scope"]["basis"])
            self.assertEqual(
                hashlib.sha256(N1_RATES.read_bytes()).hexdigest(),
                summary["pricing"]["rate_manifest_sha256"],
            )

    def test_supports_explicit_n2_rates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            evidence = self.make_evidence(root, "n2-custom-2-8192")
            output = root / "summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(evidence),
                    "--rates",
                    str(N2_RATES),
                    "--output",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))

            self.assertEqual("N2 custom", summary["pricing"]["machine_family"])
            self.assertIn("actual N2 custom", summary["cost_scope"]["basis"])
            self.assertTrue(
                all(
                    attempt["machine_type"] == "n2-custom-2-8192"
                    for attempt in summary["attempts"]
                )
            )

    def test_prices_predefined_and_custom_attempts_in_both_markets(self) -> None:
        rates = REPO_ROOT / "config/backends/gcp/gcp-rates-n1-us-west2-20260911.json"
        for predefined_attempt in (1, 2):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                evidence = self.make_evidence(root, "custom-2-12288")
                path = evidence / "batch-jobs" / f"job-test-{predefined_attempt}.json"
                job = json.loads(path.read_text())
                job["status"]["taskGroups"]["group0"]["instances"][0]["machineType"] = "n1-highmem-2"
                job["allocationPolicy"]["instances"][0]["policy"]["machineType"] = "n1-highmem-2"
                self.write_json(path, job)
                self.write_manifest(evidence)
                output = root / "summary.json"
                command = [sys.executable, str(SCRIPT), str(evidence), "--rates", str(rates), "--output", str(output)]
                subprocess.run(command, check=True, capture_output=True, text=True)
                result = json.loads(output.read_text())
                self.assertEqual(["n1-highmem-2"], result["pricing"]["predefined_machine_types"])
                for row in result["attempts"]:
                    predefined = row["attempt"] == predefined_attempt
                    ram = 13 if predefined else 12
                    market = row["provisioning_model"]
                    cpu_rate, ram_rate = {
                        (True, "SPOT"): (0.01675, 0.002243),
                        (True, "STANDARD"): (0.03797, 0.005089),
                        (False, "SPOT"): (0.01758, 0.002355),
                        (False, "STANDARD"): (0.0398685, 0.00534345),
                    }[(predefined, market)]
                    hours = row["phase_seconds"]["batch_running"] / 3600
                    self.assertEqual(2, row["provisioned_vcpu"])
                    self.assertEqual(ram, row["provisioned_memory_gib"])
                    self.assertAlmostEqual(2 * hours * cpu_rate, row["modeled_worker_cost_usd"]["vcpu"])
                    self.assertAlmostEqual(ram * hours * ram_rate, row["modeled_worker_cost_usd"]["memory"])
                # A custom-only manifest must not silently price a predefined VM.
                command[command.index(str(rates))] = str(N1_RATES)
                failed = subprocess.run(command, capture_output=True, text=True)
                self.assertNotEqual(0, failed.returncode)
                self.assertIn("machine type does not match", failed.stderr)

    def test_accepts_declared_zero_runtime_infrastructure_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            evidence = self.make_evidence(root)
            metadata_path = evidence / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            failed_attempt = metadata["calls"]["rnaseq_pipeline.star_align"][0]
            failed_attempt["executionStatus"] = "RetryableFailure"
            self.write_json(metadata_path, metadata)

            job_path = evidence / "batch-jobs" / "job-test-1.json"
            job = json.loads(job_path.read_text(encoding="utf-8"))
            job["status"]["runDuration"] = "0s"
            job["status"]["statusEvents"] = [
                {
                    "description": "Job state is set from QUEUED to SCHEDULED",
                    "eventTime": "2026-08-30T00:00:02Z",
                },
                {
                    "description": (
                        "Job state is set from SCHEDULED to FAILED. Task state is "
                        "updated from PENDING to FAILED on zones/us-west1-a/instances/1 "
                        "due to VM is recreated during task execution with exit code 50006."
                    ),
                    "eventTime": "2026-08-30T00:01:45Z",
                },
            ]
            self.write_json(job_path, job)

            self.write_manifest(evidence)

            output = root / "summary.json"
            subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence), "--output", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads(output.read_text(encoding="utf-8"))
            attempt = summary["attempts"][0]
            self.assertFalse(attempt["batch_execution_started"])
            self.assertEqual(
                "vm_recreated_before_execution", attempt["pre_execution_failure"]
            )
            self.assertIsNone(attempt["monitoring"])
            self.assertEqual(0, attempt["phase_seconds"]["batch_running"])
            self.assertEqual(0, attempt["modeled_worker_cost_usd"]["total"])
            self.assertEqual(0, summary["totals"]["modeled_worker_cost_usd"]["failed_work"])

    def test_rejects_unexplained_zero_runtime(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            evidence = self.make_evidence(Path(temporary_directory))
            job_path = evidence / "batch-jobs" / "job-test-1.json"
            job = json.loads(job_path.read_text(encoding="utf-8"))
            job["status"]["runDuration"] = "0s"
            self.write_json(job_path, job)
            self.write_manifest(evidence)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(0, result.returncode)
            self.assertIn("Batch transitions contradict runDuration", result.stderr)

    def test_rejects_missing_monitor_for_unrecognized_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            evidence = self.make_evidence(Path(temporary_directory))
            job_path = evidence / "batch-jobs" / "job-test-1.json"
            job = json.loads(job_path.read_text(encoding="utf-8"))
            job["status"]["statusEvents"][-1]["description"] = (
                "Job state is set from RUNNING to FAILED"
            )
            self.write_json(job_path, job)
            self.write_manifest(evidence)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(0, result.returncode)
            self.assertIn("missing undeclared monitoring evidence", result.stderr)

    def test_accepts_legacy_bundle_without_unavailable_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            evidence = self.make_evidence(Path(temporary_directory))
            metadata_path = evidence / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["calls"]["rnaseq_pipeline.star_align"] = [
                metadata["calls"]["rnaseq_pipeline.star_align"][1]
            ]
            self.write_json(metadata_path, metadata)

            capture_path = evidence / "capture-status.json"
            capture = json.loads(capture_path.read_text(encoding="utf-8"))
            capture["attempt_count"] = 1
            del capture["expected_unavailable_artifact_count"]
            self.write_json(capture_path, capture)
            (evidence / "batch-jobs" / "job-test-1.json").unlink()
            for stream in ("stdout", "stderr"):
                (evidence / "task-streams" / f"star_align.shard-0.attempt-1.{stream}").unlink()
            (evidence / "expected-unavailable-artifacts.tsv").unlink()
            self.write_manifest(evidence)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(evidence)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(0, result.returncode, result.stderr)

    def test_machine_family_must_match_rate_manifest(self) -> None:
        cases = (
            ("N2 with default N1 rates", "n2-custom-2-8192", None, "N1 custom"),
            ("N1 with N2 rates", "custom-2-8192", N2_RATES, "N2 custom"),
            ("unsupported N2D", "n2d-custom-2-8192", None, "N1 custom"),
            ("unpriced extended memory", "custom-2-16384-ext", None, "N1 custom"),
        )
        for name, machine_type, rates, expected in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                evidence = self.make_evidence(Path(temporary), machine_type)
                command = [sys.executable, str(SCRIPT), str(evidence)]
                if rates is not None:
                    command.extend(["--rates", str(rates)])
                result = subprocess.run(
                    command,
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(0, result.returncode)
                self.assertIn(
                    f"machine type does not match {expected} rate manifest",
                    result.stderr,
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
