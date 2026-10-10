import csv
import importlib.util
import json
import copy
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sheet = load("sample_sheet", "wdl/merge_results/sample_sheet.py")
builder = load("make_merge_inputs", "scripts/make_merge_inputs.py")
QC_FIELDS = ["sample", "reads", "pct_GC", "pct_rRNA", "pct_uniquely_mapped", "pct_multimapped", "pct_coding", "pct_utr"]


class SampleSheetTests(unittest.TestCase):
    def test_flags_failed_samples_and_sizes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "metadata.tsv").write_text("sample\texpression_mode\na\tumi_molecules\nb\tall_read\nc\tumi_molecules\n")
            with (root / "qc.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(QC_FIELDS)
                writer.writerow(["a", 30e6, 50, 1, 90, 3, 40, 40])
                writer.writerow(["b", 10e6, 50, 25, 50, 3, 30, 10])
            (root / "diagnostics").mkdir()
            for sample, strand, converged in (("a", 0.98, True), ("b", 0.72, False)):
                (root / f"diagnostics/{sample}.qc_diagnostics.json").write_text(json.dumps({
                    "sample": sample, "rsem": {"converged": converged, "iterations": 5000},
                    "picard_strand": {"correct_strand_fraction": strand},
                    "feature_counts": {"assigned_alignment_fraction": 0.8}}))
            (root / "sizes.tsv").write_text("a\t3\t2\t1.5\t0.5\nb\t1\t1\t1\t1\n")
            (root / "failed.txt").write_text("c\n")
            sheet.build(root / "metadata.tsv", root / "qc.csv", root / "diagnostics", root / "sizes.tsv",
                        root / "failed.txt", root / "sheet.tsv")
            rows = {row["sample"]: row for row in csv.DictReader((root / "sheet.tsv").read_text().splitlines(), delimiter="\t")}
            self.assertEqual(["a", "b", "c"], list(rows))
            self.assertEqual(("completed", "", "3.000"), (rows["a"]["status"], rows["a"]["qc_flags"], rows["a"]["raw_fastq_gib"]))
            self.assertEqual("low_reads_after_trim;high_rrna;low_mapped;low_exonic;low_mapped_vs_cohort;low_strand;rsem_not_converged",
                             rows["b"]["qc_flags"])
            self.assertEqual(("failed", ""), (rows["c"]["status"], rows["c"]["rsem_iterations"]))
            (root / "failed.txt").write_text("")
            with self.assertRaises(ValueError):
                sheet.build(root / "metadata.tsv", root / "qc.csv", root / "diagnostics", root / "sizes.tsv",
                            root / "failed.txt", root / "sheet2.tsv")


class MergeInputTests(unittest.TestCase):
    def test_metadata_lookup_and_recovery_errors(self):
        root = "gs://bucket/rnaseq_pipeline/run1"
        metadata = {"id": "run1", "workflowRoot": root + "/", "status": "Failed"}
        saved = subprocess.CompletedProcess([], 0, json.dumps(metadata), "")
        missing = subprocess.CompletedProcess([], 1, "", "URLs matched no objects or files")
        denied = subprocess.CompletedProcess([], 1, "", "403 Permission denied")
        with mock.patch.object(builder.subprocess, "run", return_value=saved) as run:
            self.assertEqual(metadata, builder.workflow_metadata(root))
            run.assert_called_once_with(["gcloud", "storage", "cat", root + "/metadata.json"],
                                        capture_output=True, text=True)
        wrong = subprocess.CompletedProcess([], 0, json.dumps({**metadata, "workflowRoot": "gs://other/run1"}), "")
        cases = [([denied], RuntimeError, "403 Permission denied"),
                 ([missing, denied], RuntimeError, "cannot recover metadata"),
                 ([missing, FileNotFoundError()], RuntimeError, "recovery requires caper"),
                 ([missing, wrong], ValueError, "does not match run root")]
        for responses, error, message in cases:
            with self.subTest(message=message), mock.patch.object(builder.subprocess, "run", side_effect=responses) as run:
                with self.assertRaisesRegex(error, message):
                    builder.workflow_metadata(root)
                self.assertEqual(len(responses), run.call_count)

    def test_latest_complete_results_and_failed_samples(self):
        run, retry = "gs://bucket/rnaseq_pipeline/run1", "gs://bucket/rnaseq_pipeline/run2"
        def metadata(samples, calls, release="rn8_v116"):
            document = {"status": "Failed", "inputs": {"sample_prefix": samples, "reference_release": release}, "calls": {}}
            for task, shard, status, outputs in calls:
                document["calls"].setdefault("rnaseq_pipeline." + task, []).append(
                    {"shardIndex": shard, "executionStatus": status, "outputs": outputs})
            return document
        def results(root, shard, sample, mode, rsem_status="Done"):
            rsem = "umi_molecule_rsem" if mode == "umi" else "rsem_quant"
            fc = "umi_molecule_feature_counts_task" if mode == "umi" else "feature_counts"
            path = lambda task, name: f"{root}/call-{task}/shard-{shard}/{sample}{name}"
            return [(rsem, shard, rsem_status, {"genes": path(rsem, ".genes.results"), "isoforms": path(rsem, ".isoforms.results"),
                                               "convergence": path(rsem, ".rsem_convergence.tsv")}),
                    (fc, shard, "Done", {"fc_out": path(fc, ".out")}),
                    ("qc_report", shard, "Done", {"rnaseq_report": path("qc_report", "_qc_info.csv"),
                                                 "diagnostics": path("qc_report", ".qc_diagnostics.json")}),
                    ("star_align", shard, "Done", {"bam_file": path("star_align", ".Aligned.sortedByCoord.out.bam"),
                                                   "transcriptome_bam": path("star_align", ".Aligned.toTranscriptome.out.bam")})]
        runs = {run: metadata(["a", "b", "c"], results(run, 0, "a", "umi") + results(run, 1, "b", "all", "RetryableFailure") +
                              [("udup", 0, "Done", {"molecule_transcriptome_bam": [f"{run}/call-udup/shard-0/a.umi_molecules.transcriptome.bam"]}),
                               ("qc_report", 2, "Done", {"rnaseq_report": f"{run}/call-qc_report/shard-2/c_qc_info.csv"})]),
                retry: metadata(["b"], [call for call in results(retry, 0, "b", "all") if call[0] != "feature_counts"] +
                                # cache hit: the featureCounts output stays in run1
                                [("feature_counts", 0, "Done", {"fc_out": f"{run}/call-feature_counts/shard-1/b.out"})])}
        sizes = {f"gs://raw/{s}_{m}.fastq.gz": 2**30 for s in "abc" for m in ("R1", "R2", "I1")}
        sizes.update({f"{run}/call-star_align/shard-0/a.Aligned.sortedByCoord.out.bam": 2**31,
                      f"{run}/call-star_align/shard-0/a.Aligned.toTranscriptome.out.bam": 2**30,
                      f"{run}/call-udup/shard-0/a.umi_molecules.transcriptome.bam": 2**29,
                      f"{retry}/call-star_align/shard-0/b.Aligned.sortedByCoord.out.bam": 2**30})
        inputs = {"rnaseq_pipeline." + key: value for key, value in {
            "sample_prefix": ["a", "b", "c"], "fastq1": [f"gs://raw/{s}_R1.fastq.gz" for s in "abc"],
            "fastq2": [f"gs://raw/{s}_R2.fastq.gz" for s in "abc"],
            "fastq_index": ["gs://raw/a_I1.fastq.gz", "", "gs://raw/c_I1.fastq.gz"],
            "reference_release": "rn8_v116", "output_report_name": "cohort", "merge_results_docker": "image"}.items()}
        for key in builder.SCIENTIFIC_INPUTS | {"attach_umi_docker", "umi_dup_docker"}:
            inputs["rnaseq_pipeline." + key] = "same-" + key
        for record in runs.values():
            selected = record["inputs"]["sample_prefix"]
            record["inputs"].update(builder.resolved_inputs(inputs))
            record["inputs"]["sample_prefix"] = selected
            for key in ("fastq1", "fastq2", "fastq_index"):
                record["inputs"][key] = [inputs["rnaseq_pipeline." + key]["abc".index(sample)] for sample in selected]
        runs[retry]["inputs"]["star_ramGB"] = 72  # Resource changes remain compatible.
        list_sizes = lambda *urls: {url: sizes[url] for url in urls if url in sizes}
        document = builder.build(inputs, [run, retry], runs.get, list_sizes)
        responses = []
        for root in (run, retry):
            metadata = {**runs[root], "id": root.rsplit("/", 1)[-1], "workflowRoot": root}
            responses.extend([subprocess.CompletedProcess([], 1, "", "URLs matched no objects or files"),
                              subprocess.CompletedProcess([], 0, json.dumps(metadata), "")])
        with mock.patch.object(builder.subprocess, "run", side_effect=responses) as commands:
            self.assertEqual(document, builder.build(inputs, [run, retry], list_sizes=list_sizes))
            self.assertEqual(["caper", "metadata", "run1"], commands.call_args_list[1].args[0])
            self.assertEqual(["caper", "metadata", "run2"], commands.call_args_list[3].args[0])
        empty = {root: {**record, "calls": {}} for root, record in runs.items()}
        with self.assertRaisesRegex(ValueError, "no completed samples"):
            builder.build(inputs, [run, retry], empty.get, list_sizes)
        get = lambda key: document["rnaseq_merge." + key]
        self.assertEqual((["a", "b"], ["c"]), (get("sample_prefix"), get("failed_samples")))
        self.assertEqual([["a", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"],
                          ["b", "rn8_v116", "0", "1", "all_read", "skipped_no_umi"],
                          ["c", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"]], get("expression_metadata_rows")[1:])
        self.assertEqual(f"{retry}/call-rsem_quant/shard-0/b.genes.results", get("rsem_gene_results")[1])
        self.assertEqual(f"{run}/call-feature_counts/shard-1/b.out", get("feature_counts_files")[1])
        self.assertEqual(["a", "3.000000", "2.000000", "1.000000", "0.500000"], get("sample_size_rows")[0])
        self.assertEqual(["b", "2.000000", "1.000000", "", ""], get("sample_size_rows")[1])
        self.assertEqual(builder.MERGE_RESULTS_DOCKER, get("merge_results_docker"))
        self.assertIn(["b", "fc", f"{run}/call-feature_counts/shard-1/b.out"], get("source_rows"))
        self.assertIn(f"{retry}/call-rsem_quant/shard-0/b.rsem_convergence.tsv", get("supporting_files"))
        unknown_sizes = builder.build(inputs, [run, retry], runs.get, lambda *urls: {})
        self.assertEqual(["a", "", "", "", ""], unknown_sizes["rnaseq_merge.sample_size_rows"][0])
        for key, value in (("reference_release", "rn7"), ("fastq1", ["gs://other/a", "b", "c"]),
                           ("use_umi_molecule_expression", False), ("trim_trailing_i1_base", True),
                           ("gtf_file", "other.gtf"), ("rsem_docker", "other-image")):
            changed = copy.deepcopy(runs)
            changed[run]["inputs"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                builder.build(inputs, [run, retry], changed.get, list_sizes)
        runs[retry]["status"] = "Running"
        with self.assertRaisesRegex(ValueError, "finished workflow"):
            builder.build(inputs, [run, retry], runs.get, list_sizes)


if __name__ == "__main__":
    unittest.main()
