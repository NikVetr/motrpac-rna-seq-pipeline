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
        def results(root, shard, sample, task, attempt=""):
            prefix = f"{root}/call-{{}}/shard-{shard}/{attempt}"
            rsem = "umi_molecule_rsem" if task == "umi" else "rsem_quant"
            fc = "umi_molecule_feature_counts_task" if task == "umi" else "feature_counts"
            return {prefix.format(rsem) + f"rsem_reference/{sample}.genes.results": 1,
                    prefix.format(rsem) + f"rsem_reference/{sample}.isoforms.results": 1,
                    prefix.format(fc) + f"{sample}.out": 1, prefix.format("qc_report") + f"{sample}_qc_info.csv": 1,
                    prefix.format("qc_report") + f"{sample}.qc_diagnostics.json": 1,
                    prefix.format("star_align") + f"star_out/{sample}.Aligned.sortedByCoord.out.bam": 2**31,
                    prefix.format("star_align") + f"star_out/{sample}.Aligned.toTranscriptome.out.bam": 2**30}
        listings = {run + "/": {**results(run, 0, "a", "umi"), **results(run, 1, "b", "all"),
                                f"{run}/call-qc_report/shard-2/c_qc_info.csv": 1,
                                f"{run}/call-udup/shard-0/a.umi_molecules.transcriptome.bam": 2**29},
                    retry + "/": results(retry, 1, "b", "all", "attempt-2/")}
        raw = {f"gs://raw/{s}_{m}.fastq.gz": 2**30 for s in "abc" for m in ("R1", "R2", "I1")}
        list_sizes = lambda *urls: listings[urls[0]] if urls[0] in listings else {u: raw[u] for u in urls}
        inputs = {"rnaseq_pipeline." + key: value for key, value in {
            "sample_prefix": ["a", "b", "c"], "fastq1": [f"gs://raw/{s}_R1.fastq.gz" for s in "abc"],
            "fastq2": [f"gs://raw/{s}_R2.fastq.gz" for s in "abc"],
            "fastq_index": ["gs://raw/a_I1.fastq.gz", "", "gs://raw/c_I1.fastq.gz"],
            "reference_release": "rn8_v116", "output_report_name": "cohort", "merge_results_docker": "image"}.items()}
        document = builder.build(inputs, [run, retry], list_sizes)
        get = lambda key: document["rnaseq_merge." + key]
        self.assertEqual((["a", "b"], ["c"]), (get("sample_prefix"), get("failed_samples")))
        self.assertEqual([["a", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"],
                          ["b", "rn8_v116", "0", "1", "all_read", "skipped_no_umi"],
                          ["c", "rn8_v116", "1", "0", "umi_molecules", "deduplicated"]], get("expression_metadata_rows")[1:])
        self.assertTrue(get("rsem_gene_results")[1].startswith(retry + "/call-rsem_quant/shard-1/attempt-2/"))
        self.assertEqual(["a", "3.000000", "2.000000", "1.000000", "0.500000"], get("sample_size_rows")[0])
        self.assertEqual(["b", "2.000000", "2.000000", "1.000000", "1.000000"], get("sample_size_rows")[1])


if __name__ == "__main__":
    unittest.main()
