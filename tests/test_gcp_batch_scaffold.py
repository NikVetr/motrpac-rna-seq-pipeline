import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class GcpBatchScaffoldTests(unittest.TestCase):
    def test_batch_backend_replaces_retired_papi_for_benchmarks(self) -> None:
        config = (
            REPO_ROOT / "config/backends/gcp/google_batch.conf"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory",
            config,
        )
        self.assertIn("MOTRPAC_GCP_PROJECT", config)
        self.assertIn("MOTRPAC_GCP_BATCH_ROOT", config)
        self.assertIn("MOTRPAC_GCP_BATCH_LOCATION", config)
        self.assertIn("MOTRPAC_GCP_COMPUTE_SERVICE_ACCOUNT", config)
        self.assertNotIn('compute-service-account = "default"', config)
        self.assertIn("abort-jobs-on-terminate = true", config)
        self.assertIn("max-concurrent-workflows = 1", config)
        self.assertIn("max-workflow-launch-count = 1", config)
        self.assertIn("max-scatter-width-per-scatter = 1", config)
        self.assertIn("total-max-jobs-per-root-workflow = 25", config)
        self.assertIn("concurrent-job-limit = 3", config)
        self.assertIn("maximum-polling-interval = 60", config)
        self.assertIn("batch-timeout = 4 hours", config)
        self.assertNotIn("batch-timeout = 7 days", config)
        self.assertIn("virtual-private-cloud {", config)
        self.assertIn('network-name = "default"', config)
        self.assertIn(
            'subnetwork-name = "projects/${projectId}/regions/*/subnetworks/default"',
            config,
        )
        self.assertIn("call-caching {\n  enabled = false", config)
        self.assertNotIn("PipelinesApiLifecycleActorFactory", config)
        self.assertNotIn("genomics.endpoint-url", config)
        self.assertNotIn("reference-disk-localization-manifests", config)

    def test_cromwell_release_is_pinned(self) -> None:
        manifest = json.loads(
            (
                REPO_ROOT / "config/backends/gcp/cromwell-release.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(1, manifest["schema_version"])
        self.assertEqual("92", manifest["cromwell"]["version"])
        self.assertEqual(
            "https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar",
            manifest["cromwell"]["jar_url"],
        )
        self.assertEqual(
            "e0e3a050d4124e81369a79059e5774142b2f06bd89df4a0b035f559db85cedf5",
            manifest["cromwell"]["sha256"],
        )
        self.assertEqual("Eclipse Temurin", manifest["java"]["distribution"])
        self.assertEqual(17, manifest["java"]["major_version"])

    def test_benchmark_options_are_cold_and_on_demand(self) -> None:
        options = json.loads(
            (
                REPO_ROOT
                / "config/backends/gcp/workflow-options-benchmark.example.json"
            ).read_text(encoding="utf-8")
        )
        self.assertFalse(options["read_from_cache"])
        self.assertFalse(options["write_to_cache"])
        self.assertFalse(options["use_reference_disks"])
        self.assertEqual(
            "gs://omicspipelines-public-resources/rnaseq/monitoring/"
            "sha256-9fcd1f7179d1d4106c39a68e3f116bb57f645b847158d5e1673efcf21adc2f4f/"
            "monitor_resources.sh",
            options["monitoring_script"],
        )
        self.assertNotIn("user_service_account_json", options)

    def test_benchmark_assets_are_generation_and_checksum_pinned(self) -> None:
        assets = json.loads(
            (
                REPO_ROOT / "config/backends/gcp/benchmark-assets-v1.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(1, assets["schema_version"])
        self.assertEqual("published", assets["publication_state"])
        monitor = assets["monitoring_script"]
        self.assertIn("sha256-{}".format(monitor["sha256"]), monitor["uri"])
        self.assertRegex(monitor["generation"], r"^[1-9][0-9]+$")

        sets = assets["controlled_inputs"]["sets"]
        self.assertEqual({"100k", "5m"}, set(sets))
        for expected_records, controlled in (
            (100_000, sets["100k"]),
            (5_000_000, sets["5m"]),
        ):
            self.assertEqual(expected_records, controlled["records"])
            self.assertTrue(controlled["prefix"].startswith("gs://omicspipelines-get/"))
            for artifact in [controlled["manifest"], *controlled["fastqs"].values()]:
                self.assertTrue(artifact["uri"].startswith(controlled["prefix"]))
                self.assertRegex(artifact["generation"], r"^[1-9][0-9]+$")
                self.assertRegex(artifact["sha256"], r"^[0-9a-f]{64}$")
                self.assertGreater(artifact["size_bytes"], 0)

    def test_monitoring_script_syntax_and_header(self) -> None:
        script = REPO_ROOT / "scripts/gcp/monitor_resources.sh"
        subprocess.run(["bash", "-n", str(script)], check=True)
        environment = os.environ.copy()
        environment["MOTRPAC_MONITOR_INTERVAL_SECONDS"] = "0"
        environment["MOTRPAC_MONITOR_MAX_SAMPLES"] = "1"
        result = subprocess.run(
            ["bash", str(script)],
            check=True,
            capture_output=True,
            encoding="utf-8",
            env=environment,
        )
        lines = result.stdout.splitlines()
        self.assertEqual(2, len(lines))
        self.assertEqual(9, len(lines[0].split("\t")))
        self.assertEqual(9, len(lines[1].split("\t")))
        self.assertTrue(
            re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", lines[1].split("\t")[0])
        )

    def test_existing_infrastructure_preflight_is_read_only(self) -> None:
        script = REPO_ROOT / "scripts/gcp/preflight_existing_infrastructure.sh"
        subprocess.run(["bash", "-n", str(script)], check=True)
        contents = script.read_text(encoding="utf-8")
        self.assertNotIn("set -e", contents)
        self.assertIn("run_check", contents)
        self.assertIn("READ_ONLY_PREFLIGHT_FAILURES", contents)
        self.assertIn('--location="$batch_location"', contents)
        for mutating_command in (
            "batch jobs submit",
            "builds submit",
            "compute instances start",
            "storage cp",
        ):
            self.assertNotIn(mutating_command, contents)

    def test_running_vm_watcher_is_read_only(self) -> None:
        script = REPO_ROOT / "scripts/gcp/watch_running_vms.sh"
        subprocess.run(["bash", "-n", str(script)], check=True)
        contents = script.read_text(encoding="utf-8")
        self.assertIn("labels", contents)
        self.assertIn("batch-node", contents)
        self.assertIn("batch-job-id", contents)
        self.assertIn("goog-batch-worker", contents)
        for mutating_command in (
            "batch jobs submit",
            "compute instances create",
            "compute instances delete",
            "compute instances start",
            "compute instances stop",
        ):
            self.assertNotIn(mutating_command, contents)

    def test_evidence_capture_preserves_every_attempt(self) -> None:
        script = REPO_ROOT / "scripts/gcp/capture_workflow_evidence.sh"
        subprocess.run(["bash", "-n", str(script)], check=True)

        workflow_id = "12345678-1234-1234-1234-123456789abc"
        revision = subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], text=True
        ).strip()
        with tempfile.TemporaryDirectory() as temporary_directory:
            temp = Path(temporary_directory)
            fake_bin = temp / "bin"
            fake_bin.mkdir()
            fake_gcloud = fake_bin / "gcloud"
            fake_gcloud.write_text(
                """#!/bin/sh
set -eu
if [ "$1 $2 $3" = "batch jobs describe" ]; then
    printf '{"name":"%s"}\\n' "$4"
elif [ "$1 $2 $3" = "storage objects describe" ]; then
    object_name=${4#gs://test/}
    printf '%s\\n' '{"bucket":"test","name":"'"$object_name"'","generation":"123","metageneration":"1","size":"42","md5_hash":"YWJj","crc32c_hash":"ZGVm"}'
elif [ "$1 $2" = "storage cp" ]; then
    printf 'captured %s\\n' "$3" >"$4"
else
    exit 2
fi
""",
                encoding="utf-8",
            )
            fake_gcloud.chmod(0o755)

            attempts = []
            for attempt in (1, 2):
                stem = "job-test-{}".format(attempt)
                attempts.append(
                    {
                        "attempt": attempt,
                        "executionStatus": "Done" if attempt == 2 else "Failed",
                        "jobId": "projects/test-project/locations/us-west1/jobs/{}".format(
                            stem
                        ),
                        "monitoringLog": "gs://test/{}/monitoring.log".format(stem),
                        "preemptible": attempt == 1,
                        "shardIndex": 0,
                        "stderr": "gs://test/{}/stderr".format(stem),
                        "stdout": "gs://test/{}/stdout".format(stem),
                    }
                )
            metadata = temp / "metadata.json"
            metadata.write_text(
                json.dumps(
                    {
                        "id": workflow_id,
                        "status": "Succeeded",
                        "submittedFiles": {
                            "inputs": json.dumps(
                                {
                                    "fastq": "gs://test/input_R1.fastq.gz",
                                    "reference": "gs://test/reference.tar.gz",
                                }
                            ),
                            "options": json.dumps({"option": "value"}),
                        },
                        "outputs": {
                            "rnaseq_pipeline.rsem_genes_count": "gs://test/rsem_genes_count.txt",
                            "rnaseq_pipeline.umi_metrics": [
                                "gs://test/sample.umi_metrics.json"
                            ],
                        },
                        "calls": {"rnaseq_pipeline.star_align": attempts},
                    }
                ),
                encoding="utf-8",
            )
            output = temp / "evidence"
            environment = os.environ.copy()
            environment["PATH"] = "{}:{}".format(fake_bin, environment["PATH"])
            subprocess.run(
                ["bash", str(script), str(metadata), str(output), revision],
                check=True,
                capture_output=True,
                encoding="utf-8",
                env=environment,
            )

            status = json.loads(
                (output / "capture-status.json").read_text(encoding="utf-8")
            )
            self.assertEqual(workflow_id, status["workflow_id"])
            self.assertEqual(2, status["attempt_count"])
            self.assertEqual(2, status["submitted_gcs_object_count"])
            self.assertEqual(2, status["top_level_output_object_count"])
            self.assertEqual(0, status["missing_artifact_count"])
            self.assertTrue(status["complete"])
            repository = json.loads(
                (output / "repository.json").read_text(encoding="utf-8")
            )
            self.assertRegex(repository["revision"], r"^[0-9a-f]{40}$")
            self.assertEqual(
                repository["revision"], repository["expected_submission_revision"]
            )
            self.assertIsInstance(repository["clean"], bool)
            objects = json.loads(
                (output / "input-objects.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                ["gs://test/input_R1.fastq.gz", "gs://test/reference.tar.gz"],
                [entry["uri"] for entry in objects],
            )
            self.assertEqual({42}, {int(entry["size_bytes"]) for entry in objects})
            output_objects = json.loads(
                (output / "output-objects.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                {
                    "rnaseq_pipeline.rsem_genes_count",
                    "rnaseq_pipeline.umi_metrics",
                },
                {entry["output_name"] for entry in output_objects},
            )
            self.assertEqual(
                2, len(list((output / "top-level-outputs").iterdir()))
            )
            self.assertEqual(2, len(list((output / "batch-jobs").glob("*.json"))))
            self.assertEqual(6, len(list((output / "task-streams").iterdir())))
            self.assertTrue((output / "evidence-manifest.sha256").is_file())

            mismatch = subprocess.run(
                ["bash", str(script), str(metadata), str(temp / "wrong"), "b" * 40],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertNotEqual(0, mismatch.returncode)
            self.assertIn("does not match the submitted revision", mismatch.stderr)

    def test_operator_docs_match_capture_and_runtime_limits(self) -> None:
        runbook = (REPO_ROOT / "docs/gcp-cli-canary-runbook.md").read_text(
            encoding="utf-8"
        )
        normalized_runbook = " ".join(runbook.split())
        self.assertIn("at most three concurrent Batch workers", normalized_runbook)
        self.assertIn("scripts/gcp/capture_workflow_evidence.sh", runbook)
        self.assertIn("local user account", normalized_runbook)
        self.assertIn(
            "controller's attached `cromwell-prod` identity", normalized_runbook
        )


if __name__ == "__main__":
    unittest.main()
