import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MERGE_ROOT = REPO_ROOT / "wdl/merge_results"


class MergeResultsTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.fc = self.root / "fc"
        self.rsem = self.root / "rsem"
        self.qc = self.root / "qc"
        for directory in (self.fc, self.rsem, self.qc):
            directory.mkdir()
        self.samples = ["sample.beta", "sample.alpha.1"]
        self.order = self.root / "sample_order.txt"
        self.order.write_text("\n".join(self.samples) + "\n", encoding="utf-8")

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write_fixtures(self):
        for sample, values in {
            "sample.alpha.1": ("11", "12"),
            "sample.beta": ("21", "22"),
        }.items():
            (self.fc / (sample + ".out")).write_text(
                "# Program:featureCounts\n"
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\tinput.bam\n"
                "gene.two\tchr1\t1\t2\t+\t2\t{}\n"
                "gene.one\tchr1\t3\t4\t+\t2\t{}\n".format(*values),
                encoding="utf-8",
            )
            with (self.rsem / (sample + ".genes.results")).open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "gene_id",
                        "transcript_id(s)",
                        "length",
                        "effective_length",
                        "expected_count",
                        "TPM",
                        "FPKM",
                    ]
                )
                writer.writerow(["gene.two", "tx2", "2", "1", values[0], "1.2", "2.3"])
                writer.writerow(["gene.one", "tx1", "2", "1", values[1], "3.4", "4.5"])
            (self.qc / (sample + "_qc_info.csv")).write_text(
                "sample,reads,pct_GC,pct_umi_dup,pct_chimeric\n"
                "{},{},52.5,,\n".format(sample, values[0]),
                encoding="utf-8",
            )

    def run_script(self, script, *arguments):
        return subprocess.run(
            [sys.executable, str(MERGE_ROOT / script)] + list(arguments),
            cwd=self.root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_shuffled_files_preserve_declared_samples_dotted_ids_and_gene_order(self):
        self.write_fixtures()
        commands = [
            (
                "merge_fc.py",
                "--fc-dir",
                str(self.fc),
                "--sample-order",
                str(self.order),
            ),
            (
                "merge_rsem.py",
                "--rsem-dir",
                str(self.rsem),
                "--sample-order",
                str(self.order),
            ),
            (
                "consolidate_qc_report.py",
                "--qc-dir",
                str(self.qc),
                "--sample-order",
                str(self.order),
                "--output-name",
                "cohort.csv",
            ),
        ]
        for command in commands:
            result = self.run_script(*command)
            self.assertEqual(0, result.returncode, result.stderr.decode())

        feature_counts = (self.root / "featureCounts.txt").read_text().splitlines()
        self.assertEqual("gene_id\tsample.beta\tsample.alpha.1", feature_counts[0])
        self.assertEqual("gene.two\t21\t11", feature_counts[1])
        rsem_counts = (self.root / "rsem_genes_count.txt").read_text().splitlines()
        self.assertEqual("gene_id\tsample.beta\tsample.alpha.1", rsem_counts[0])
        self.assertEqual("gene.two\t21\t11", rsem_counts[1])
        with (self.root / "cohort.csv").open(encoding="utf-8", newline="") as handle:
            qc_rows = list(csv.DictReader(handle))
        self.assertEqual(self.samples, [row["sample"] for row in qc_rows])
        self.assertTrue(all(row["pct_chimeric"] == "" for row in qc_rows))

    def test_sample_set_mismatch_fails_loudly(self):
        self.write_fixtures()
        (self.fc / "sample.alpha.1.out").unlink()
        result = self.run_script(
            "merge_fc.py",
            "--fc-dir",
            str(self.fc),
            "--sample-order",
            str(self.order),
        )
        self.assertNotEqual(0, result.returncode)

    def test_optional_qc_blanks_merge_but_core_blanks_fail(self):
        self.write_fixtures()
        for sample in self.samples:
            path = self.qc / (sample + "_qc_info.csv")
            path.write_text(
                path.read_text(encoding="utf-8").replace(",52.5,,\n", ",,,\n"),
                encoding="utf-8",
            )
        result = self.run_script(
            "consolidate_qc_report.py",
            "--qc-dir",
            str(self.qc),
            "--sample-order",
            str(self.order),
            "--output-name",
            "optional_blanks.csv",
        )
        self.assertEqual(0, result.returncode, result.stderr.decode())

        path = self.qc / "sample.beta_qc_info.csv"
        path.write_text(
            path.read_text(encoding="utf-8").replace("sample.beta,21,", "sample.beta,,"),
            encoding="utf-8",
        )
        failed = self.run_script(
            "consolidate_qc_report.py",
            "--qc-dir",
            str(self.qc),
            "--sample-order",
            str(self.order),
            "--output-name",
            "core_blank.csv",
        )
        self.assertNotEqual(0, failed.returncode)


if __name__ == "__main__":
    unittest.main()
