import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/gcp"))
import capture_task_profiles as profiles
from summarize_workflow_cost import monitor, verify_evidence


class TaskProfileTests(unittest.TestCase):
    def test_sample_and_release_labels_accept_resolved_and_submitted_inputs(self):
        for qualified, submitted in ((False, False), (True, False), (True, True)):
            with self.subTest(qualified=qualified, submitted=submitted), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                prefix = "rnaseq_pipeline." if qualified else ""
                inputs = {prefix + "sample_prefix": ["rat_a", "rat_b"], prefix + "reference_release": "rn8_v116"}
                metadata = {"id": "workflow", "calls": {"rnaseq_pipeline.star_align": [
                    {"shardIndex": 1, "attempt": 1, "executionStatus": "Done", "inputs": {"ncpu": 12}}]}}
                if submitted:
                    metadata["submittedFiles"] = {"inputs": json.dumps(inputs)}
                else:
                    metadata["inputs"] = inputs
                path = root / "metadata.json"
                path.write_text(json.dumps(metadata))
                row, = profiles.collect(path, root / "capture")["attempts"]
                self.assertEqual("rat_b", row["sample"])
                self.assertEqual("rn8_v116", row["reference_release"])

    def test_working_memory_uses_simultaneous_samples_and_accepts_old_logs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "monitoring"
            header = "epoch_s\tcpu_usage_usec\tmemory_peak_bytes\tdisk_used_kb\tdisk_available_kb\tmemory_current_bytes"
            path.write_text(header + "\tmemory_inactive_file_bytes\tmemory_anon_bytes\n"
                            "1\t0\tNA\t0\t1\t10737418240\t9663676416\t1073741824\n"
                            "2\t1\tNA\t0\t1\t8589934592\t1073741824\t6442450944\n")
            result = monitor(path)
            self.assertEqual(7, result["peak_working_set_gib"])
            self.assertEqual(6, result["peak_memory_anon_gib"])
            path.write_text(header + "\n1\t0\tNA\t0\t1\t1073741824\n")
            self.assertIsNone(monitor(path)["peak_working_set_gib"])

    def test_failed_workflow_retains_attempts_counts_and_large_object_sizes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.json"
            metrics = json.dumps({"denominators": {"selected_representative_qnames_present_in_transcriptome": 7}}).encode()
            objects = {"gs://b/input.bam": None, "gs://b/counts.isoforms.results": None,
                       "gs://b/s.rsem_convergence.tsv": b"component\ttranscript_id\n" + b"transcript\tt1\n" * 100000,
                       "gs://b/s.rsem_gene_convergence.tsv": b"gene_id\titeration\n",
                       "gs://b/s.umi_molecule_expression_metrics.json": metrics}
            metadata.write_text(json.dumps({"id": "workflow", "status": "Failed", "calls": {
                "rnaseq_pipeline.udup": [
                    {"shardIndex": 0, "attempt": 1, "executionStatus": "Failed", "jobId": "job1",
                     "inputs": {"star_align": "gs://b/input.bam"}, "failures": [{"message": "OOM"}]},
                    {"shardIndex": 1, "attempt": 1, "executionStatus": "Done", "jobId": "job2",
                     "inputs": {"star_align": "gs://b/input.bam"},
                     "outputs": {"results": "gs://b/counts.isoforms.results",
                                 "convergence": "gs://b/s.rsem_convergence.tsv",
                                 "gene_convergence": "gs://b/s.rsem_gene_convergence.tsv",
                                 "molecule_expression_metrics": ["gs://b/s.umi_molecule_expression_metrics.json"]}},
                    {"shardIndex": 2, "attempt": 1, "executionStatus": "Running"}]} }))
            downloads = []
            def cloud(*args):
                if args[:3] == ("storage", "objects", "describe"):
                    data = objects[args[3]]
                    return json.dumps({"size": len(data) if data else 10**12, "generation": "123"}).encode()
                if args[:3] == ("batch", "jobs", "describe"):
                    return json.dumps({"name": args[3], "status": {"state": "FAILED" if args[3] == "job1" else "SUCCEEDED"}}).encode()
                self.assertEqual(("storage", "cat"), args[:2])
                downloads.append(args[2])
                self.assertTrue(args[2].endswith("#123"))
                data = objects[args[2].split("#")[0]]
                self.assertIsNotNone(data, "large data must not be downloaded")
                return data
            with patch.object(profiles, "cloud", side_effect=cloud):
                result = profiles.collect(metadata, root / "capture")
            self.assertEqual(2, len(result["attempts"]))
            self.assertEqual(1, len(result["not_terminal"]))
            self.assertEqual(10**12, result["attempts"][0]["input_bytes"]["star_align"])
            self.assertEqual("OOM", result["attempts"][0]["failures"][0]["message"])
            self.assertEqual(3, len(downloads))
            self.assertEqual(objects["gs://b/s.rsem_convergence.tsv"].decode(),
                             result["attempts"][1]["metrics"]["s.rsem_convergence.tsv"])
            self.assertEqual("gene_id\titeration\n", result["attempts"][1]["metrics"]["s.rsem_gene_convergence.tsv"])
            self.assertEqual(7, result["attempts"][1]["metrics"][
                "s.umi_molecule_expression_metrics.json"]["denominators"][
                "selected_representative_qnames_present_in_transcriptome"])
            self.assertTrue(result["complete_capture"])
            verify_evidence(root / "capture")
            (root / "capture/profiles.json").write_text("tampered")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                verify_evidence(root / "capture")

    def test_unavailable_input_is_retained_as_incomplete_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.json"
            metadata.write_text(json.dumps({"id": "workflow", "calls": {"rnaseq_pipeline.rsem_quant": [
                {"shardIndex": 0, "attempt": 1, "executionStatus": "Failed",
                 "inputs": {"transcriptome_bam": "gs://b/missing.bam"}}]}}))
            with patch.object(profiles, "cloud", side_effect=subprocess.CalledProcessError(1, "gcloud")):
                with self.assertRaisesRegex(ValueError, "incomplete capture"):
                    profiles.collect(metadata, root / "capture")
            result = json.loads((root / "capture/profiles.json").read_text())
            self.assertFalse(result["complete_capture"])
            self.assertEqual(1, len(result["attempts"]))
            self.assertEqual("gs://b/missing.bam", result["errors"][0]["uri"])
            verify_evidence(root / "capture")

    def test_only_preemption_404_is_nonfatal_and_retry_logs_are_not_substituted(self):
        for status, error, fatal in (("RetryableFailure", "404 not found", False),
                                     ("Done", "404 not found", True),
                                     ("RetryableFailure", "403 permission denied for shard-404/monitoring.log", True)):
            with self.subTest(status=status, error=error), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                path = root / "metadata.json"
                path.write_text(json.dumps({"id": "workflow", "calls": {"rnaseq_pipeline.udup": [
                    {"shardIndex": 0, "attempt": 1, "executionStatus": status,
                     "failures": [{"message": "VMPreemption(50001)"}], "monitoringLog": "gs://b/monitoring.log"},
                    {"shardIndex": 0, "attempt": 2, "executionStatus": "Done",
                     "monitoringLog": "gs://b/attempt-2/monitoring.log"}]}}))
                trace = b"epoch_s\tcpu_usage_usec\tmemory_peak_bytes\tdisk_used_kb\tdisk_available_kb\tmemory_current_bytes\n1\t0\t1024\t0\t1\t1024\n"
                def cloud(*args):
                    if args[:3] == ("storage", "objects", "describe"):
                        if args[3] == "gs://b/monitoring.log":
                            raise ValueError(error)
                        return json.dumps({"size": len(trace), "generation": "123"}).encode()
                    self.assertEqual(("storage", "cat", "gs://b/attempt-2/monitoring.log#123"), args)
                    return trace
                with patch.object(profiles, "cloud", side_effect=cloud):
                    if fatal:
                        with self.assertRaisesRegex(ValueError, "incomplete capture"):
                            profiles.collect(path, root / "capture")
                    else:
                        profiles.collect(path, root / "capture")
                result = json.loads((root / "capture/profiles.json").read_text())
                self.assertFalse(result["complete_capture"])
                self.assertIsNone(result["attempts"][0]["monitoring"])
                self.assertIsNotNone(result["attempts"][1]["monitoring"])
                self.assertEqual(not fatal, result["errors"][0]["expected_after_preemption"])
                verify_evidence(root / "capture")

    def test_batch_permission_error_still_preserves_other_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "metadata.json"
            path.write_text(json.dumps({"id": "workflow", "calls": {"rnaseq_pipeline.udup": [
                {"shardIndex": 0, "attempt": 1, "executionStatus": "Done", "jobId": "job1",
                 "inputs": {"ncpu": 4}}]}}))
            with patch.object(profiles, "cloud", side_effect=ValueError("PERMISSION_DENIED")):
                with self.assertRaisesRegex(ValueError, "incomplete capture"):
                    profiles.collect(path, root / "capture")
            result = json.loads((root / "capture/profiles.json").read_text())
            self.assertEqual("job1", result["errors"][0]["job"])
            self.assertEqual(1, len(result["attempts"]))
            verify_evidence(root / "capture")


if __name__ == "__main__":
    unittest.main()
