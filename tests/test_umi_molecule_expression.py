from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
UMI_DIR = REPO_ROOT / "wdl" / "umi_dup"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


propagate = load_module(
    "propagate_molecule_qnames", UMI_DIR / "propagate_molecule_qnames.py"
)
summarize = load_module(
    "summarize_molecule_expression", UMI_DIR / "summarize_molecule_expression.py"
)


class FakeRecord:
    def __init__(self, query_name: str, rx: str | None):
        self.query_name = query_name
        self.rx = rx

    def has_tag(self, name: str) -> bool:
        return name == "RX" and self.rx is not None

    def get_tag(self, name: str) -> str:
        assert name == "RX" and self.rx is not None
        return self.rx


class UmiMoleculeExpressionTests(unittest.TestCase):
    def test_representative_umi_contract_is_strict(self) -> None:
        record = FakeRecord("read:AAAAAAAA", "acgtacgt")
        self.assertEqual("ACGTACGT", propagate.umi_from_record(record, "rx_v1", 8))
        for invalid in (None, "ACGT", "ACGTNCGT", "ACGTXCGT"):
            with self.assertRaises(ValueError):
                propagate.umi_from_record(
                    FakeRecord("read:AAAAAAAA", invalid), "rx_v1", 8
                )
        with self.assertRaisesRegex(ValueError, "unsupported UMI representation"):
            propagate.umi_from_record(record, "legacy_qname", 8)

    def test_projection_denominators_fail_loudly(self) -> None:
        qnames = {
            "selected_representative_qnames": 5,
            "selected_genomic_present_qnames": 5,
            "selected_genomic_absent_qnames": 0,
            "selected_transcriptome_present_qnames": 4,
            "selected_transcriptome_absent_qnames": 1,
            "transcriptome_qnames_absent_from_genomic_source": 0,
            "selected_genomic_alignment_records": 13,
            "selected_transcriptome_alignment_records": 14,
        }
        genomic = {"kept_alignment_records": 13}
        transcriptome = {"kept_alignment_records": 14}
        propagate.validate_denominators(qnames, genomic, transcriptome)

        qnames["selected_transcriptome_absent_qnames"] = 0
        with self.assertRaisesRegex(ValueError, "present/absent counts disagree"):
            propagate.validate_denominators(qnames, genomic, transcriptome)

    def test_summary_reconciles_grouping_and_projection_denominators(self) -> None:
        umi = {
            "algorithm": "umi_tools_directional",
            "acceptance_class": "scientific-truth",
            "umi_molecules": 5,
            "umi_eligible_acgt_templates": 10,
            "tool_version": "UMI-tools version: 1.1.6",
            "container": "image@sha256:" + "a" * 64,
            "edit_distance": 1,
            "random_seed": 12345,
            "representation": "rx_v1",
        }
        propagation = {
            "status": "production-primary",
            "algorithm": "propagate_genomic_umi_representative_qnames_v1",
            "runtime_versions": {"pysam": "0.22.1", "sqlite": "3.47.2"},
            "qnames": {
                "selected_representative_qnames": 5,
                "selected_genomic_alignment_records": 13,
                "transcriptome_source_unique_qnames": 9,
                "selected_transcriptome_present_qnames": 4,
                "selected_transcriptome_absent_qnames": 1,
                "selected_transcriptome_alignment_records": 14,
            },
            "genomic": {"kept_alignment_records": 13},
            "transcriptome": {"kept_alignment_records": 14},
        }
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            umi_path = directory / "umi.json"
            propagation_path = directory / "propagation.json"
            output_path = directory / "output.json"
            umi_path.write_text(json.dumps(umi), encoding="utf-8")
            propagation_path.write_text(json.dumps(propagation), encoding="utf-8")
            arguments = ["summarize_molecule_expression.py"]
            arguments += ["--umi-metrics", str(umi_path)]
            arguments += ["--propagation-metrics", str(propagation_path)]
            arguments += ["--output", str(output_path)]
            with patch.object(sys, "argv", arguments):
                self.assertEqual(0, summarize.main())
            result = json.loads(output_path.read_text(encoding="utf-8"))
            propagation["qnames"]["selected_representative_qnames"] = 4
            propagation_path.write_text(json.dumps(propagation), encoding="utf-8")
            with patch.object(sys, "argv", arguments):
                with self.assertRaisesRegex(ValueError, "counts differ"):
                    summarize.main()

        self.assertEqual(
            "umi_tools_directional_molecule_expression_v1",
            result["algorithm"],
        )
        self.assertEqual("production-primary", result["status"])
        self.assertEqual(
            "directional UMI molecule expression",
            result["interpretation"]["canonical_outputs"],
        )
        denominators = result["denominators"]
        self.assertEqual(
            4, denominators["selected_representative_qnames_present_in_transcriptome"]
        )
        self.assertEqual(
            1, denominators["selected_representative_qnames_absent_from_transcriptome"]
        )

    def test_workflow_makes_molecule_expression_primary_and_all_read_optional(
        self,
    ) -> None:
        workflow = (REPO_ROOT / "wdl/rnaseq_pipeline_scatter.wdl").read_text()
        umi_wdl = (UMI_DIR / "umi_dup.wdl").read_text()
        dockerfile = (REPO_ROOT / "dockerfiles/umi_dup.Dockerfile").read_text()
        merge = (REPO_ROOT / "wdl/merge_results/merge_expression.wdl").read_text()

        for expected in (
            "--no-sort-output",
            "--extract-umi-method=tag",
            "--umi-tag=RX",
            "--method=directional",
            "--edit-distance-threshold=1",
            "--multimapping-detection-method=NH",
            "--random-seed=12345",
            ">/dev/null",
            "Boolean emit_molecule_expression = false",
            "Array[File] molecule_genomic_bam",
            "Array[File] molecule_transcriptome_bam",
        ):
            self.assertIn(expected, umi_wdl)
        self.assertEqual(1, umi_wdl.count("--method=directional"))
        self.assertNotIn("nudup", umi_wdl.lower())
        self.assertIn("File umi_metrics", umi_wdl)
        self.assertIn("propagate_molecule_qnames.py", umi_wdl)
        self.assertIn("summarize_molecule_expression.py", umi_wdl)
        self.assertEqual(
            "FROM quay.io/biocontainers/umi_tools@sha256:"
            "94c7cd9a713157affe93d3f1fa60e60d35a6385adc6b419d5f73c68eea8a54e8",
            dockerfile.splitlines()[0],
        )
        for script in (
            "prepare_umi_bam.py",
            "summarize_umi_tools.py",
            "propagate_molecule_qnames.py",
            "summarize_molecule_expression.py",
        ):
            self.assertIn("COPY wdl/umi_dup/{}".format(script), dockerfile)
        self.assertNotIn("python:2", dockerfile)
        for retired in ("nudup.py", "umi_dup.sh"):
            self.assertFalse((UMI_DIR / retired).exists())

        for expected in (
            "Boolean use_umi_molecule_expression = true",
            "Boolean retain_all_read_expression = false",
            "Boolean run_all_read_expression =",
            "if (run_all_read_expression)",
            "umi_report=udup.umi_report",
            "Array[File] umi_metrics = select_all(udup.umi_metrics)",
            "if !use_umi_molecule_expression || has_fastq_index then [true] else []",
            "Boolean umi_expression_inputs_valid = umi_expression_input_contract[0]",
            "transcriptome_align=if use_umi_molecule_expression then "
            "[star_align.transcriptome_bam] else []",
            "rsem_files=primary_rsem_genes",
            "feature_counts_files=primary_feature_counts",
            "rsem_report=primary_rsem_report",
            "fc_report=primary_feature_counts_report",
            "if (use_umi_molecule_expression && retain_all_read_expression)",
            'output_prefix="all_read"',
            "rsem_files=select_all(rsem_quant.genes)",
            "feature_counts_files=select_all(feature_counts.fc_out)",
        ):
            self.assertIn(expected, workflow)
        molecule_calls = workflow[
            workflow.index("call fc.feature_counts as umi_molecule_feature_counts_task") :
            workflow.index("File primary_rsem_genes")
        ]
        self.assertEqual(2, molecule_calls.count("SID=sample_prefix[i]"))
        self.assertNotIn('SID=sample_prefix[i] + ".umi_molecules"', molecule_calls)
        output_block = workflow.rsplit("output {", 1)[1]
        self.assertNotIn("molecule_genomic_bam", output_block)
        self.assertNotIn("molecule_transcriptome_bam", output_block)
        self.assertIn("Array[File] umi_expression_metrics", output_block)
        for output in (
            "all_read_rsem_genes_count",
            "all_read_rsem_genes_tpm",
            "all_read_rsem_genes_fpkm",
            "all_read_feature_counts",
        ):
            self.assertIn(f"File? {output}", output_block)
        self.assertNotIn("File? umi_molecule_rsem", output_block)
        self.assertNotIn("File? umi_molecule_feature", output_block)

        self.assertIn("/usr/local/src/merge_rsem.py", merge)
        self.assertIn("/usr/local/src/merge_fc.py", merge)
        self.assertIn("String output_prefix", merge)
        self.assertNotIn("call expression_merge.merge_expression as merge_umi_expression", workflow)
        self.assertNotIn("consolidate_qc_report.py", merge)


if __name__ == "__main__":
    unittest.main()
