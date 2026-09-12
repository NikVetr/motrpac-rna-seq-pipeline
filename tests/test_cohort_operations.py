import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


locality = load("locality", "scripts/gcp/check_locality.py")
metadata = load("metadata", "scripts/prepare_sample_metadata.py")


class CohortOperationsTests(unittest.TestCase):
    def test_shutdown_requires_complete_checksummed_local_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = root / "evidence"
            evidence.mkdir()
            status = evidence / "capture-status.json"
            status.write_text(json.dumps({"complete": False, "workflow_status": "Succeeded",
                                          "missing_artifact_count": 0}))
            gcloud = root / "gcloud"
            gcloud.write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$CALL_LOG"\n'
                              'if [ "$3" = describe ]; then echo TERMINATED; fi\n')
            gcloud.chmod(0o755)
            log = root / "calls"
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ["PATH"], CALL_LOG=str(log))
            command = ["bash", str(ROOT / "scripts/gcp/stop_controller_after_capture.sh"),
                       str(evidence), "project", "controller", "us-west2-a"]
            self.assertNotEqual(0, subprocess.run(command, env=env, capture_output=True).returncode)
            self.assertFalse(log.exists())
            status.write_text(json.dumps({"complete": True, "workflow_status": "Succeeded",
                                          "missing_artifact_count": 0}))
            # Missing checksum evidence must also prevent a stop.
            self.assertNotEqual(0, subprocess.run(command, env=env, capture_output=True).returncode)
            self.assertFalse(log.exists())
            manifest = subprocess.check_output(["sha256sum", "capture-status.json"], cwd=evidence)
            (evidence / "evidence-manifest.sha256").write_bytes(manifest)
            result = subprocess.run(command, env=env, capture_output=True)
            self.assertEqual(0, result.returncode, result.stderr)
            self.assertIn("instances stop controller", log.read_text())

    def test_locality_requires_regional_execution_but_allows_remote_sources(self):
        options = {"default_runtime_attributes": {"zones": "us-west2-a us-west2-b"}}
        def describe(bucket):
            return {"location": "US-WEST2" if bucket == "execution" else "US",
                    "locationType": "region" if bucket == "execution" else "multi-region"}
        report = locality.check("gs://execution/run", "us-west2", options,
                                {"fastqs": ["gs://raw/sample.fastq.gz"]}, describe)
        self.assertEqual([True, False], [row["colocated"] for row in report])
        with self.assertRaisesRegex(ValueError, "execution bucket"):
            locality.check("gs://raw/run", "us-west2", options, {}, describe)
        with self.assertRaisesRegex(ValueError, "zones differ"):
            locality.check("gs://execution/run", "us-west1", options, {}, describe)
        options["gcp_batch_gcs_root"] = "gs://elsewhere/run"
        with self.assertRaisesRegex(ValueError, "execution root"):
            locality.check("gs://execution/run", "us-west2", options, {}, describe)

    def test_metadata_preserves_order_values_and_rejects_missing_or_old_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            matrix, study, qc, output = [root / name for name in ("matrix", "study", "qc", "output")]
            matrix.write_text("gene_id\tb\ta\ngene\t1\t2\n")
            study.write_text("sample,pid,RIN\na,001,9.40\nb,002,8.0\nc,003,7\n")
            qc.write_text("sample,pct_umi_dup\na,80\nb,90\n")
            metadata.prepare(matrix, study, qc, output, ["pid", "RIN"])
            self.assertEqual("sample,pid,RIN,pct_umi_dup\nb,002,8.0,90\na,001,9.40,80\n", output.read_text())
            with self.assertRaises(FileExistsError):
                metadata.prepare(matrix, study, qc, output, ["pid"])
            study.write_text("sample,pid,RIN\na,001,NA\nb,002,8\n")
            with self.assertRaisesRegex(ValueError, "missing required"):
                metadata.prepare(matrix, study, qc, root / "missing", ["RIN"])
            study.write_text("sample,pid,pct_umi_dup\na,001,75\nb,002,85\n")
            with self.assertRaisesRegex(ValueError, "overlap"):
                metadata.prepare(matrix, study, qc, root / "old", ["pid"])
            study.write_text("sample,pid\na,001\na,002\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                metadata.prepare(matrix, study, qc, root / "duplicate", ["pid"])

    def test_predefined_policy_rejects_other_regions_and_platform_overrides(self):
        inputs = {"rnaseq_pipeline.prefer_predefined_n1": True}
        def describe(bucket):
            return {"location": "US-WEST2", "locationType": "region"}
        options = {"default_runtime_attributes": {"zones": "us-west2-a"}}
        locality.check("gs://execution/run", "us-west2", options, inputs, describe)
        for key, value in (("cpuPlatform", "Intel Cascade Lake"), ("predefinedMachineType", "e2-highmem-2")):
            conflicting = {"default_runtime_attributes": dict(options["default_runtime_attributes"], **{key: value})}
            with self.assertRaisesRegex(ValueError, "predefined N1 policy"):
                locality.check("gs://execution/run", "us-west2", conflicting, inputs, describe)
        with self.assertRaisesRegex(ValueError, "predefined N1 policy"):
            locality.check("gs://execution/run", "us-west1",
                           {"default_runtime_attributes": {"zones": "us-west1-a"}}, inputs, describe)

    def test_e2_policy_rejects_conflicting_machine_selection(self):
        inputs = {"rnaseq_pipeline.use_e2": True}
        options = {"default_runtime_attributes": {"zones": "us-west2-a"}}
        def describe(bucket):
            return {"location": "US-WEST2", "locationType": "region"}
        locality.check("gs://execution/run", "us-west2", options, inputs, describe)
        for key, value in (("cpuPlatform", "Intel Ice Lake"), ("predefinedMachineType", "e2-highmem-2")):
            conflicting = {"default_runtime_attributes": dict(options["default_runtime_attributes"], **{key: value})}
            with self.assertRaisesRegex(ValueError, "E2 policy"):
                locality.check("gs://execution/run", "us-west2", conflicting, inputs, describe)
        with self.assertRaisesRegex(ValueError, "E2 policy"):
            locality.check("gs://execution/run", "us-west2", options,
                           dict(inputs, **{"rnaseq_pipeline.prefer_predefined_n1": True}), describe)


if __name__ == "__main__":
    unittest.main()
