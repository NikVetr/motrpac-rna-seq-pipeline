from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import make_json_rnaseq as generator  # noqa: E402


class FakeGcsFileSystem:
    def __init__(self, paths: list[str], existing: set[str]) -> None:
        self.paths = paths
        self.existing = existing
        self.project = None

    def glob(self, pattern: str) -> list[str]:
        self.pattern = pattern
        return list(self.paths)

    def exists(self, path: str) -> bool:
        return path in self.existing


class InputJsonGeneratorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def arguments(self, *, include_index: bool = False) -> argparse.Namespace:
        return argparse.Namespace(
            gcp_path="gs://example/fastq_raw",
            output_path=str(self.temp),
            output_report_name="cohort.v1.csv",
            undetermined=False,
            organism="human",
            version="gencode_v39",
            num_chunks=1,
            docker_repo="registry.example/rnaseq/",
            index=include_index,
            umi_molecule_expression=False,
            project="local-test-only",
        )

    @staticmethod
    def fake_gcsfs(filesystem: FakeGcsFileSystem) -> types.ModuleType:
        module = types.ModuleType("gcsfs")

        def constructor(*, project=None):
            filesystem.project = project
            return filesystem

        module.GCSFileSystem = constructor
        return module

    def test_filtering_precedes_nonempty_deterministic_batching(self) -> None:
        batches = generator.build_batches(
            [
                "bucket/Undetermined_lane_R1.fastq.gz",
                "bucket/sample_b_R1.fastq.gz",
                "bucket/sample_a_R1.fastq.gz",
            ],
            2,
        )
        self.assertEqual(
            [["sample_a"], ["sample_b"]],
            [batch["sample_prefix"] for batch in batches],
        )
        self.assertTrue(all(batch["i1"] is None for batch in batches))
        with self.assertRaisesRegex(ValueError, "no R1 FASTQs"):
            generator.build_batches(["bucket/Undetermined_R1.fastq.gz"], 1)
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            generator.build_batches(["bucket/sample_R1.fastq.gz"], 2)

    def test_sample_and_uri_identity_are_unique(self) -> None:
        batch = generator.build_batches(
            ["bucket/sample_R1.fastq.gz"], 1, include_index=True
        )[0]
        self.assertEqual(["gs://bucket/sample_I1.fastq.gz"], batch["i1"])
        with self.assertRaisesRegex(ValueError, "duplicate sample prefix"):
            generator.build_batches(
                ["bucket/a/sample_R1.fastq.gz", "bucket/b/sample_R1.fastq.gz"], 1
            )
        with self.assertRaisesRegex(ValueError, "FASTQ URIs must be unique"):
            generator.make_json_dict(
                "human",
                "gencode_v39",
                "registry.example/rnaseq",
                "cohort",
                ["gs://example/sample.fastq.gz"],
                ["gs://example/sample.fastq.gz"],
                None,
                ["sample"],
            )
        with self.assertRaisesRegex(ValueError, "not filename-safe"):
            generator.build_batches(["bucket/not safe_R1.fastq.gz"], 1)

    def test_document_requires_aligned_nonempty_arrays(self) -> None:
        document = generator.make_json_dict(
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort.csv",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
        )
        self.assertIsNone(document["rnaseq_pipeline.fastq_index"])
        self.assertEqual("cohort", document["rnaseq_pipeline.output_report_name"])
        with self.assertRaisesRegex(ValueError, "arrays must be nonempty and aligned"):
            generator.make_json_dict(
                "human", "gencode_v39", "registry.example/rnaseq", "cohort"
            )

    def test_umi_molecule_expression_is_explicit_and_requires_i1(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires a matched I1"):
            generator.make_json_dict(
                "human",
                "gencode_v39",
                "registry.example/rnaseq",
                "cohort",
                ["gs://example/sample_R1.fastq.gz"],
                ["gs://example/sample_R2.fastq.gz"],
                None,
                ["sample"],
                use_umi_molecule_expression=True,
            )

        document = generator.make_json_dict(
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            ["gs://example/sample_I1.fastq.gz"],
            ["sample"],
            use_umi_molecule_expression=True,
        )
        self.assertTrue(document["rnaseq_pipeline.use_umi_molecule_expression"])

        default_document = generator.make_json_dict(
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
        )
        self.assertNotIn(
            "rnaseq_pipeline.use_umi_molecule_expression", default_document
        )

    def test_qc_toggles_are_default_on_and_only_disabled_values_are_emitted(self) -> None:
        arguments = (
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
        )
        default_document = generator.make_json_dict(*arguments)
        toggle_keys = {
            f"rnaseq_pipeline.{name}"
            for name in (
                "run_pretrim_fastqc",
                "run_posttrim_fastqc",
                "run_contamination_qc",
                "run_alignment_qc",
                "run_umi_qc",
            )
        }
        self.assertTrue(toggle_keys.isdisjoint(default_document))

        selective_document = generator.make_json_dict(
            *arguments,
            run_pretrim_fastqc=False,
            run_contamination_qc=False,
            run_umi_qc=False,
        )
        self.assertFalse(selective_document["rnaseq_pipeline.run_pretrim_fastqc"])
        self.assertFalse(selective_document["rnaseq_pipeline.run_contamination_qc"])
        self.assertFalse(selective_document["rnaseq_pipeline.run_umi_qc"])
        self.assertNotIn("rnaseq_pipeline.run_posttrim_fastqc", selective_document)
        self.assertNotIn("rnaseq_pipeline.run_alignment_qc", selective_document)

    def test_contamination_qc_fusion_and_sampling_are_opt_in(self) -> None:
        arguments = (
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
        )
        default_document = generator.make_json_dict(*arguments)
        self.assertNotIn(
            "rnaseq_pipeline.combine_contamination_qc", default_document
        )
        self.assertNotIn("rnaseq_pipeline.contamination_qc_pairs", default_document)

        combined_document = generator.make_json_dict(
            *arguments, combine_contamination_qc=True
        )
        self.assertTrue(
            combined_document["rnaseq_pipeline.combine_contamination_qc"]
        )
        self.assertNotIn("rnaseq_pipeline.contamination_qc_pairs", combined_document)

        sampled_document = generator.make_json_dict(
            *arguments, contamination_qc_pairs=1_000_000
        )
        self.assertEqual(
            1_000_000,
            sampled_document["rnaseq_pipeline.contamination_qc_pairs"],
        )
        self.assertNotIn(
            "rnaseq_pipeline.combine_contamination_qc", sampled_document
        )

        with self.assertRaisesRegex(ValueError, "nonnegative integer"):
            generator.make_json_dict(*arguments, contamination_qc_pairs=-1)
        with self.assertRaisesRegex(ValueError, "requires run_contamination_qc"):
            generator.make_json_dict(
                *arguments,
                run_contamination_qc=False,
                combine_contamination_qc=True,
            )

    def test_main_checks_mates_before_writing(self) -> None:
        r1 = "bucket/sample_R1.fastq.gz"
        r2 = "gs://bucket/sample_R2.fastq.gz"
        filesystem = FakeGcsFileSystem([r1], {r2})
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(filesystem)}):
            self.assertEqual(0, generator.main(self.arguments()))
        document = json.loads((self.temp / "set1_rnaseq.json").read_text())
        self.assertEqual(["sample"], document["rnaseq_pipeline.sample_prefix"])
        self.assertEqual("local-test-only", filesystem.project)

        (self.temp / "set1_rnaseq.json").unlink()
        missing = FakeGcsFileSystem([r1], set())
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(missing)}):
            with self.assertRaisesRegex(ValueError, "required mate/index objects are missing"):
                generator.main(self.arguments())
        self.assertFalse((self.temp / "set1_rnaseq.json").exists())

    def test_existing_generated_output_fails_before_listing(self) -> None:
        (self.temp / "set2_rnaseq.json").write_text("{}", encoding="utf-8")
        filesystem = FakeGcsFileSystem([], set())
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(filesystem)}):
            with self.assertRaisesRegex(ValueError, "already contains"):
                generator.main(self.arguments())
        self.assertFalse(hasattr(filesystem, "pattern"))

    def test_organism_and_reference_are_required_together(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "scripts" / "make_json_rnaseq.py"),
                "-g",
                "gs://example/fastq_raw",
                "-o",
                str(self.temp),
                "-r",
                "cohort",
                "-n",
                "1",
                "-a",
                "human",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(2, result.returncode)
        self.assertIn(b"--version", result.stderr)


if __name__ == "__main__":
    unittest.main()
