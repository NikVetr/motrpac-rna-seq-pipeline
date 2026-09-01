from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
UMI_DIR = REPO_ROOT / "wdl" / "umi_dup"
SUMMARIZE = UMI_DIR / "summarize_umi_tools.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


prepare = load_module("prepare_umi_bam", UMI_DIR / "prepare_umi_bam.py")


class UmiGroupingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_inputs(self, *, eligible: int = 9, mates_never_found: int = 0) -> None:
        policy = {
            "schema_version": 1,
            "representation": "rx_v1",
            "umi_length": 8,
            "mapping_eligibility": "all_primary_proper",
            "primary_proper_templates": 10,
            "excluded_n_umi_templates": 1,
            "eligible_acgt_templates": eligible,
            "eligible_nh1_templates": 7,
            "eligible_multimapped_templates": 2,
        }
        (self.temp / "policy.json").write_text(json.dumps(policy))
        (self.temp / "umi.log").write_text(
            "\n".join(
                (
                    "Reads: Input Reads: 9, Read pairs: 9",
                    "Number of reads out: 6",
                    "Total number of positions deduplicated: 5",
                    "Mean number of unique UMIs per position: 1.20",
                    "Max. number of unique UMIs per position: 3",
                    "Searching for mates for 0 unmatched alignments",
                    f"{mates_never_found} mates never found",
                )
            )
            + "\n"
        )

    def summarize(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(SUMMARIZE),
                "--sample",
                "sample.1",
                "--policy",
                "policy.json",
                "--log",
                "umi.log",
                "--metrics",
                "metrics.json",
                "--report",
                "report.tsv",
                "--version",
                "UMI-tools version: 1.1.6",
                "--container",
                "example.invalid/umi_tools@sha256:abc",
            ],
            cwd=self.temp,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_qname_umi_is_strict_and_case_normalized(self) -> None:
        self.assertEqual("ACGTACGT", prepare.umi_from_query_name("read:acgtacgt"))
        self.assertEqual("ACGTNCGT", prepare.umi_from_query_name("read:acgtncgt"))
        for name in ("read", "read:ACGT", "read:ACGTXCGT", "read:"):
            with self.assertRaises(ValueError):
                prepare.umi_from_query_name(name)

    def test_summary_emits_denominator_bound_metrics_and_qc_shape(self) -> None:
        self.write_inputs()
        result = self.summarize()
        self.assertEqual(0, result.returncode, result.stderr)

        metrics = json.loads((self.temp / "metrics.json").read_text())
        self.assertEqual("umi_tools_directional", metrics["algorithm"])
        self.assertEqual("all_primary_proper", metrics["mapping_eligibility"])
        self.assertEqual(9, metrics["umi_eligible_acgt_templates"])
        self.assertEqual(6, metrics["umi_molecules"])
        self.assertEqual(3, metrics["umi_duplicate_templates"])
        self.assertAlmostEqual(100 / 3, metrics["pct_umi_dup_eligible"])
        self.assertEqual(
            "Sample\tpct_umi_dup\nsample.1\t33.333333\n",
            (self.temp / "report.tsv").read_text(),
        )

    def test_summary_rejects_denominator_or_mate_failure(self) -> None:
        self.write_inputs(eligible=10)
        result = self.summarize()
        self.assertNotEqual(0, result.returncode)
        self.assertIn("eligible policy denominator", result.stderr)

        self.write_inputs(mates_never_found=1)
        result = self.summarize()
        self.assertNotEqual(0, result.returncode)
        self.assertIn("could not find 1 mates", result.stderr)

    def test_wdl_and_image_pin_one_directional_policy(self) -> None:
        wdl = (UMI_DIR / "umi_dup.wdl").read_text()
        dockerfile = (REPO_ROOT / "dockerfiles/umi_dup.Dockerfile").read_text()
        workflow = (REPO_ROOT / "wdl/rnaseq_pipeline_scatter.wdl").read_text()

        for expected in (
            "--no-sort-output",
            "--extract-umi-method=tag",
            "--umi-tag=RX",
            "--method=directional",
            "--edit-distance-threshold=1",
            "--multimapping-detection-method=NH",
            "--random-seed=12345",
            ">/dev/null",
        ):
            self.assertIn(expected, wdl)
        self.assertEqual(1, wdl.count("--method=directional"))
        self.assertNotIn("nudup", wdl.lower())
        self.assertIn("File umi_metrics", wdl)
        self.assertIn("Boolean emit_molecule_expression = false", wdl)
        self.assertIn("propagate_molecule_qnames.py", wdl)
        self.assertIn("summarize_molecule_expression.py", wdl)
        self.assertIn("Array[File] molecule_genomic_bam", wdl)
        self.assertIn("Array[File] molecule_transcriptome_bam", wdl)

        self.assertEqual(
            "FROM quay.io/biocontainers/umi_tools@sha256:"
            "94c7cd9a713157affe93d3f1fa60e60d35a6385adc6b419d5f73c68eea8a54e8",
            dockerfile.splitlines()[0],
        )
        self.assertIn("COPY wdl/umi_dup/prepare_umi_bam.py", dockerfile)
        self.assertIn("COPY wdl/umi_dup/summarize_umi_tools.py", dockerfile)
        self.assertIn("COPY wdl/umi_dup/propagate_molecule_qnames.py", dockerfile)
        self.assertIn(
            "COPY wdl/umi_dup/summarize_molecule_expression.py", dockerfile
        )
        self.assertNotIn("python:2", dockerfile)
        self.assertFalse((UMI_DIR / "nudup.py").exists())
        self.assertFalse((UMI_DIR / "umi_dup.sh").exists())

        self.assertIn("umi_report=udup.umi_report", workflow)
        self.assertIn("Array[File] umi_metrics = select_all(udup.umi_metrics)", workflow)
        self.assertIn("Boolean use_umi_molecule_expression = true", workflow)
        self.assertIn(
            "call fc.feature_counts as umi_molecule_feature_counts_task", workflow
        )
        self.assertIn("call rsem.rsem as umi_molecule_rsem", workflow)


if __name__ == "__main__":
    unittest.main()
