import csv
import io
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
COLLECTOR = REPO_ROOT / "wdl/collect_qc_metrics/rnaseq_qc.py"


class NativeQcTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.sample = "sample.1"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def fastqc_archive(self, name, reads, filenames):
        archive_path = self.root / name
        with tarfile.open(archive_path, "w:gz") as archive:
            for index, filename in enumerate(filenames, start=1):
                data = (
                    "Filename\t{}\n"
                    "Total Sequences\t{}\n"
                    "%GC\t{}\n"
                    "#Total Deduplicated Percentage\t{}\n"
                ).format(filename, reads, 50 + index, 80 - index).encode()
                payload = io.BytesIO()
                with zipfile.ZipFile(payload, "w") as zipped:
                    zipped.writestr("mate{}/fastqc_data.txt".format(index), data)
                member = tarfile.TarInfo("reports/mate{}_fastqc.zip".format(index))
                member.size = len(payload.getvalue())
                archive.addfile(member, io.BytesIO(payload.getvalue()))
        return archive_path

    def write(self, name, content):
        path = self.root / name
        path.write_text(content, encoding="utf-8")
        return path

    def fixtures(self):
        raw_names = ("raw_R1.fastq.gz", "raw_R2.fastq.gz")
        trim_names = ("trim_R1.fastq.gz", "trim_R2.fastq.gz")
        paths = {
            "pre": self.fastqc_archive("pre.tar.gz", 100, raw_names),
            "post": self.fastqc_archive("post.tar.gz", 90, trim_names),
            "cutadapt": self.write(
                "cutadapt.log",
                "Total read pairs processed: 100\n"
                "  Read 1 with adapter: 20 (20.0%)\n"
                "  Read 2 with adapter: 40 (40.0%)\n"
                "Pairs written (passing filters): 90 (90.0%)\n"
                "Total basepairs processed: 1,000 bp\n"
                "Total written (filtered): 900 bp (90.0%)\n",
            ),
            "mapped": self.write(
                "mapped.txt",
                "Sample\tpct_chrX\tpct_chrY\tpct_chrM\tpct_chrAuto\tpct_contig\n"
                "{}\t2.12345\t0.123456\t3.4567\t94.1\t0.2\n".format(self.sample),
            ),
            "rrna": self.write("rrna.txt", "Sample\tpct_rRNA\n{}\t0.10%\n".format(self.sample)),
            "globin": self.write("globin.txt", "Sample\tpct_globin\n{}\t1.20%\n".format(self.sample)),
            "phix": self.write("phix.txt", "Sample\tpct_phix\n{}\t0.00%\n".format(self.sample)),
            "umi": self.write("umi.txt", "Sample\tpct_umi_dup\n{}\t45.12\n".format(self.sample)),
        }
        star_values = {
            "Number of input reads": "90",
            "Average input read length": "150",
            "Uniquely mapped reads number": "80",
            "Uniquely mapped reads %": "88.89%",
            "Average mapped length": "149.5",
            "Number of splices: Total": "20",
            "Number of splices: Annotated (sjdb)": "19",
            "Number of splices: GT/AG": "18",
            "Number of splices: GC/AG": "1",
            "Number of splices: AT/AC": "0",
            "Number of splices: Non-canonical": "1",
            "% of reads mapped to multiple loci": "5.00%",
            "% of reads mapped to too many loci": "1.00%",
            "% of reads unmapped: too many mismatches": "0.00%",
            "% of reads unmapped: too short": "4.00%",
            "% of reads unmapped: other": "1.11%",
            "% of chimeric reads": "0.00%",
        }
        paths["star"] = self.write(
            "Log.final.out",
            "".join("{} |\t{}\n".format(key, value) for key, value in star_values.items()),
        )
        paths["markdup"] = self.write(
            "markdup.txt",
            "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
            "PERCENT_DUPLICATION\tREAD_PAIRS_EXAMINED\n0.12345\t80\n\n",
        )
        paths["rna"] = self.write(
            "rna.txt",
            "## METRICS CLASS\tpicard.analysis.RnaSeqMetrics\n"
            "PCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\t"
            "PCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tMEDIAN_5PRIME_TO_3PRIME_BIAS\n"
            "0.5\t0.2\t0.1\t0.2\t0.7\t0.25\n\n",
        )
        return raw_names, trim_names, paths

    def run_collector(
        self, output, include_umi=True, expected_raw_names=None, omitted=()
    ):
        raw_names, trim_names, paths = self.fixtures()
        if expected_raw_names is not None:
            raw_names = expected_raw_names
        command = [
            sys.executable,
            str(COLLECTOR),
            "--sample", self.sample,
            "--pretrim-r1-filename", raw_names[0],
            "--pretrim-r2-filename", raw_names[1],
            "--posttrim-r1-filename", trim_names[0],
            "--posttrim-r2-filename", trim_names[1],
            "--cutadapt-report", str(paths["cutadapt"]),
            "--star-log", str(paths["star"]),
            "--output", str(output),
        ]
        optional_reports = {
            "pre": ("--fastqc-pretrim", paths["pre"]),
            "post": ("--fastqc-posttrim", paths["post"]),
            "mapped": ("--mapped-report", paths["mapped"]),
            "rrna": ("--rrna-report", paths["rrna"]),
            "globin": ("--globin-report", paths["globin"]),
            "phix": ("--phix-report", paths["phix"]),
            "markdup": ("--markduplicates-metrics", paths["markdup"]),
            "rna": ("--rnaseq-metrics", paths["rna"]),
            "umi": ("--umi-report", paths["umi"]),
        }
        for name, (option, path) in optional_reports.items():
            if name not in omitted and (name != "umi" or include_umi):
                command.extend([option, str(path)])
        return subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_native_reports_emit_fixed_published_columns(self):
        output = self.root / "qc.csv"
        result = self.run_collector(output)
        self.assertEqual(0, result.returncode, result.stderr.decode())
        with output.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        self.assertEqual(40, len(reader.fieldnames))
        self.assertEqual("30", rows[0]["pct_adapter_detected"])
        self.assertEqual("10", rows[0]["pct_trimmed_bases"])
        self.assertEqual("12.345", rows[0]["pct_picard_dup"])
        self.assertEqual("45.12", rows[0]["pct_umi_dup"])
        self.assertEqual("", rows[0]["pct_chimeric"])

    def test_absent_umi_is_explicit_and_filename_mismatch_fails(self):
        output = self.root / "without_umi.csv"
        result = self.run_collector(output, include_umi=False)
        self.assertEqual(0, result.returncode, result.stderr.decode())
        with output.open(encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertIn("pct_umi_dup", row)
        self.assertEqual("", row["pct_umi_dup"])

        failed = self.run_collector(
            self.root / "bad.csv",
            expected_raw_names=("wrong_R1.fastq.gz", "raw_R2.fastq.gz"),
        )
        self.assertNotEqual(0, failed.returncode)

    def test_skipped_qc_groups_leave_only_declared_optional_fields_empty(self):
        output = self.root / "core_only.csv"
        result = self.run_collector(
            output,
            include_umi=False,
            omitted={"pre", "post", "mapped", "rrna", "globin", "phix", "markdup", "rna"},
        )
        self.assertEqual(0, result.returncode, result.stderr.decode())
        with output.open(encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual("100", row["reads_raw"])
        self.assertEqual("90", row["reads"])
        self.assertEqual("88.89", row["pct_uniquely_mapped"])
        for name in (
            "pct_GC",
            "pct_dup_sequence",
            "pct_rRNA",
            "pct_globin",
            "pct_phix",
            "pct_picard_dup",
            "pct_umi_dup",
            "pct_chrX",
            "pct_chrY",
            "pct_chrM",
            "pct_chrAuto",
            "pct_contig",
            "pct_coding",
            "pct_utr",
            "pct_intronic",
            "pct_intergenic",
            "pct_mrna",
            "median_5_3_bias",
        ):
            self.assertEqual("", row[name], name)


if __name__ == "__main__":
    unittest.main()
