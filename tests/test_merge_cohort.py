import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


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
    def test_latest_complete_results_and_failed_samples(self):
        run, retry = "gs://bucket/rnaseq_pipeline/run1", "gs://bucket/rnaseq_pipeline/run2"
        def metadata(samples, calls, release="rn8_v116"):
            document = {"inputs": {"sample_prefix": samples, "reference_release": release}, "calls": {}}
            for task, shard, status, files in calls:
                document["calls"].setdefault("rnaseq_pipeline." + task, []).append(
                    {"shardIndex": shard, "executionStatus": status, "outputs": {"files": files}})
            return document
        def results(root, shard, sample, mode, rsem_status="Done"):
            rsem = "umi_molecule_rsem" if mode == "umi" else "rsem_quant"
            fc = "umi_molecule_feature_counts_task" if mode == "umi" else "feature_counts"
            path = lambda task, name: f"{root}/call-{task}/shard-{shard}/{sample}{name}"
            return [(rsem, shard, rsem_status, [path(rsem, ".genes.results"), path(rsem, ".isoforms.results")]),
                    (fc, shard, "Done", [path(fc, ".out")]),
                    ("qc_report", shard, "Done", [path("qc_report", "_qc_info.csv"), path("qc_report", ".qc_diagnostics.json")]),
                    ("star_align", shard, "Done", [path("star_align", ".Aligned.sortedByCoord.out.bam"),
                                                   path("star_align", ".Aligned.toTranscriptome.out.bam")])]
        runs = {run: metadata(["a", "b", "c"], results(run, 0, "a", "umi") + results(run, 1, "b", "all", "RetryableFailure") +
                              [("udup", 0, "Done", [[f"{run}/call-udup/shard-0/a.umi_molecules.transcriptome.bam"]]),
                               ("qc_report", 2, "Done", [f"{run}/call-qc_report/shard-2/c_qc_info.csv"])]),
                retry: metadata(["b"], [call for call in results(retry, 0, "b", "all") if call[0] != "feature_counts"] +
                                # cache hit: the featureCounts output stays in run1
                                [("feature_counts", 0, "Done", [f"{run}/call-feature_counts/shard-1/b.out"])])}
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
        list_sizes = lambda *urls: {url: sizes[url] for url in urls if url in sizes}
        document = builder.build(inputs, [run, retry], runs.get, list_sizes)
        get = lambda key: document["rnaseq_merge." + key]
        self.assertEqual((["a", "b"], ["c"]), (get("sample_prefix"), get("failed_samples")))
        self.assertEqual([["a", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"],
                          ["b", "rn8_v116", "0", "1", "all_read", "skipped_no_umi"],
                          ["c", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"]], get("expression_metadata_rows")[1:])
        self.assertEqual(f"{retry}/call-rsem_quant/shard-0/b.genes.results", get("rsem_gene_results")[1])
        self.assertEqual(f"{run}/call-feature_counts/shard-1/b.out", get("feature_counts_files")[1])
        self.assertEqual(["a", "3.000000", "2.000000", "1.000000", "0.500000"], get("sample_size_rows")[0])
        self.assertEqual(["b", "2.000000", "1.000000", "", ""], get("sample_size_rows")[1])
        runs[retry]["inputs"]["reference_release"] = "rn7"
        with self.assertRaises(ValueError):
            builder.build(inputs, [run, retry], runs.get, list_sizes)


if __name__ == "__main__":
    unittest.main()
