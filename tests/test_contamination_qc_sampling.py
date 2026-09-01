from __future__ import annotations

import gzip
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TASK_PATH = REPO_ROOT / "wdl" / "contamination_qc" / "contamination_qc.wdl"


class ContaminationQcSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        task = TASK_PATH.read_text(encoding="utf-8")
        match = re.search(r"<<'PYTHON'\n(.*?)\nPYTHON", task, flags=re.DOTALL)
        if match is None:
            raise AssertionError("could not locate embedded paired-sampling program")
        cls.program = match.group(1)

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_fastq(self, name: str, read_names: list[str], mate: int) -> Path:
        path = self.temp / name
        with gzip.open(path, "wt", encoding="ascii") as handle:
            for index, read_name in enumerate(read_names):
                sequence = "ACGT" if index % 2 == 0 else "TGCA"
                handle.write(
                    "@{}/{}\n{}\n+\nIIII\n".format(read_name, mate, sequence)
                )
        return path

    def run_sampler(
        self,
        r1: Path,
        r2: Path,
        pairs: int,
        stem: str,
        reported_pairs: int | None = None,
    ) -> subprocess.CompletedProcess[str]:
        if reported_pairs is None:
            with gzip.open(r1, "rt", encoding="ascii") as handle:
                reported_pairs = sum(1 for _ in handle) // 4
        report = self.temp / "{}_cutadapt.log".format(stem)
        report.write_text(
            "Pairs written (passing filters): {:,} (100.0%)\n".format(
                reported_pairs
            ),
            encoding="utf-8",
        )
        return subprocess.run(
            [
                sys.executable,
                "-",
                str(r1),
                str(r2),
                str(report),
                str(pairs),
                "motrpac-contamination-qc-v1",
                str(self.temp / "{}_R1.fastq.gz".format(stem)),
                str(self.temp / "{}_R2.fastq.gz".format(stem)),
                str(self.temp / "{}_manifest.json".format(stem)),
            ],
            input=self.program,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_sample_is_exact_paired_deterministic_and_auditable(self) -> None:
        names = ["read{:02d}".format(index) for index in range(20)]
        r1 = self.write_fastq("input_R1.fastq.gz", names, 1)
        r2 = self.write_fastq("input_R2.fastq.gz", names, 2)

        first = self.run_sampler(r1, r2, 7, "first")
        second = self.run_sampler(r1, r2, 7, "second")
        self.assertEqual("", first.stderr)
        self.assertEqual(0, first.returncode)
        self.assertEqual(0, second.returncode, second.stderr)

        def sampled_names(path: Path) -> list[str]:
            with gzip.open(path, "rt", encoding="ascii") as handle:
                lines = handle.readlines()
            return [lines[index].split("/", 1)[0][1:] for index in range(0, len(lines), 4)]

        first_r1 = sampled_names(self.temp / "first_R1.fastq.gz")
        first_r2 = sampled_names(self.temp / "first_R2.fastq.gz")
        second_r1 = sampled_names(self.temp / "second_R1.fastq.gz")
        self.assertEqual(7, len(first_r1))
        self.assertEqual(first_r1, first_r2)
        self.assertEqual(first_r1, second_r1)
        self.assertEqual(first_r1, sorted(first_r1))

        manifest = json.loads(
            (self.temp / "first_manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual("sha256-counter-floyd-ordinal-v1", manifest["algorithm"])
        self.assertEqual(20, manifest["input_pairs"])
        self.assertEqual(7, manifest["selected_pairs"])
        self.assertEqual("motrpac-contamination-qc-v1", manifest["seed"])
        self.assertFalse(manifest["used_full_input"])
        self.assertRegex(manifest["selected_name_sha256"], r"^[0-9a-f]{64}$")

    def test_unsynchronized_fastqs_fail_loudly(self) -> None:
        r1 = self.write_fastq("input_R1.fastq.gz", ["a", "b"], 1)
        r2 = self.write_fastq("input_R2.fastq.gz", ["a", "wrong"], 2)
        result = self.run_sampler(r1, r2, 1, "bad")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("paired FASTQs are unsynchronized at record 2", result.stderr)

    def test_oversized_samples_use_all_pairs_and_invalid_inputs_fail(self) -> None:
        r1 = self.write_fastq("input_R1.fastq.gz", ["a", "b"], 1)
        r2 = self.write_fastq("input_R2.fastq.gz", ["a", "b"], 2)
        oversized = self.run_sampler(r1, r2, 3, "oversized")
        self.assertEqual(0, oversized.returncode, oversized.stderr)
        manifest = json.loads(
            (self.temp / "oversized_manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(3, manifest["requested_pairs"])
        self.assertEqual(2, manifest["selected_pairs"])
        self.assertTrue(manifest["used_full_input"])

        negative = self.run_sampler(r1, r2, -1, "negative")
        self.assertNotEqual(0, negative.returncode)
        self.assertIn("must be nonnegative", negative.stderr)

        count_mismatch = self.run_sampler(
            r1, r2, 1, "count_mismatch", reported_pairs=3
        )
        self.assertNotEqual(0, count_mismatch.returncode)
        self.assertIn("Cutadapt reports 3 pairs, but the FASTQs contain 2", count_mismatch.stderr)

    def test_zero_pairs_records_full_depth_without_writing_a_sample(self) -> None:
        r1 = self.write_fastq("input_R1.fastq.gz", ["a", "b"], 1)
        r2 = self.write_fastq("input_R2.fastq.gz", ["a", "b"], 2)
        result = self.run_sampler(r1, r2, 0, "full")
        self.assertEqual(0, result.returncode, result.stderr)
        self.assertFalse((self.temp / "full_R1.fastq.gz").exists())
        self.assertFalse((self.temp / "full_R2.fastq.gz").exists())
        manifest = json.loads(
            (self.temp / "full_manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual("none-full-depth", manifest["algorithm"])
        self.assertEqual(2, manifest["input_pairs"])
        self.assertEqual(2, manifest["selected_pairs"])
        self.assertTrue(manifest["used_full_input"])

    def test_wrong_illumina_read_number_fails_loudly(self) -> None:
        r1 = self.temp / "input_R1.fastq.gz"
        r2 = self.temp / "input_R2.fastq.gz"
        with gzip.open(r1, "wt", encoding="ascii") as handle:
            handle.write("@read 2:N:0:INDEX\nACGT\n+\nIIII\n")
        with gzip.open(r2, "wt", encoding="ascii") as handle:
            handle.write("@read 2:N:0:INDEX\nACGT\n+\nIIII\n")
        result = self.run_sampler(r1, r2, 1, "wrong_role")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("wrong read-number field", result.stderr)


if __name__ == "__main__":
    unittest.main()
