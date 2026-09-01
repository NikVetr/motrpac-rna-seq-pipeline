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
    "summarize_molecule_expression",
    UMI_DIR / "summarize_molecule_expression.py",
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
            "status": "production-shadow",
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
            arguments = [
                "summarize_molecule_expression.py",
                "--umi-metrics",
                str(umi_path),
                "--propagation-metrics",
                str(propagation_path),
                "--output",
                str(output_path),
            ]
            with patch.object(sys, "argv", arguments):
                self.assertEqual(0, summarize.main())
            result = json.loads(output_path.read_text(encoding="utf-8"))
            propagation["qnames"]["selected_representative_qnames"] = 4
            propagation_path.write_text(json.dumps(propagation), encoding="utf-8")
            with patch.object(sys, "argv", arguments):
                with self.assertRaisesRegex(ValueError, "counts differ"):
                    summarize.main()

        self.assertEqual(
            "umi_tools_directional_molecule_expression_shadow_v1",
            result["algorithm"],
        )
        self.assertEqual("production-shadow", result["status"])
        self.assertEqual(
            4,
            result["denominators"][
                "selected_representative_qnames_present_in_transcriptome"
            ],
        )
        self.assertEqual(
            1,
            result["denominators"][
                "selected_representative_qnames_absent_from_transcriptome"
            ],
        )

    def test_workflow_keeps_all_read_outputs_and_namespaces_molecule_matrices(
        self,
    ) -> None:
        workflow = (REPO_ROOT / "wdl/rnaseq_pipeline_scatter.wdl").read_text()
        merge = (REPO_ROOT / "wdl/merge_results/merge_expression.wdl").read_text()

        self.assertIn("Boolean use_umi_molecule_expression = false", workflow)
        self.assertIn(
            "if !use_umi_molecule_expression || has_fastq_index then [true] else []",
            workflow,
        )
        self.assertIn(
            "Boolean umi_expression_inputs_valid = umi_expression_input_contract[0]",
            workflow,
        )
        self.assertIn(
            "transcriptome_align=if use_umi_molecule_expression then "
            "[star_align.transcriptome_bam] else []",
            workflow,
        )
        self.assertIn("rsem_files=rsem_quant.genes", workflow)
        self.assertIn("feature_counts_files=feature_counts.fc_out", workflow)
        self.assertIn(
            "rsem_files=select_all(umi_molecule_rsem.genes)", workflow
        )
        self.assertIn(
            "feature_counts_files=select_all(umi_molecule_feature_counts_task.fc_out)",
            workflow,
        )
        molecule_calls = workflow[
            workflow.index("call fc.feature_counts as umi_molecule_feature_counts_task") :
            workflow.index("call collect_qc.rnaseqQC as qc_report")
        ]
        self.assertEqual(2, molecule_calls.count("SID=sample_prefix[i]"))
        self.assertNotIn('SID=sample_prefix[i] + ".umi_molecules"', molecule_calls)
        output_block = workflow.rsplit("output {", 1)[1]
        self.assertNotIn("molecule_genomic_bam", output_block)
        self.assertNotIn("molecule_transcriptome_bam", output_block)
        for output in (
            "umi_molecule_rsem_genes_count",
            "umi_molecule_rsem_genes_tpm",
            "umi_molecule_rsem_genes_fpkm",
            "umi_molecule_feature_counts",
        ):
            self.assertIn(f"File? {output}", output_block)

        self.assertIn("/usr/local/src/merge_rsem.py", merge)
        self.assertIn("/usr/local/src/merge_fc.py", merge)
        self.assertNotIn("consolidate_qc_report.py", merge)


if __name__ == "__main__":
    unittest.main()
