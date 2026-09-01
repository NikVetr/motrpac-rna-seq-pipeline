from __future__ import annotations

import argparse
import json
import re
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

    def arguments(self, *, include_index: bool = True) -> argparse.Namespace:
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
            umi_molecule_expression=include_index,
            star_disk_type=None,
            project="local-test-only",
        )

    @staticmethod
    def document_arguments(*, include_index: bool = False) -> tuple:
        return (
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            ["gs://example/sample_I1.fastq.gz"] if include_index else None,
            ["sample"],
        )

    def make_document(self, *, include_index: bool = False, **options) -> dict:
        return generator.make_json_dict(
            *self.document_arguments(include_index=include_index), **options
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
                *self.document_arguments()[:4],
                ["gs://example/sample.fastq.gz"],
                ["gs://example/sample.fastq.gz"],
                None,
                ["sample"],
            )
        with self.assertRaisesRegex(ValueError, "not filename-safe"):
            generator.build_batches(["bucket/not safe_R1.fastq.gz"], 1)

    def test_document_requires_aligned_nonempty_arrays(self) -> None:
        arguments = list(self.document_arguments())
        arguments[3] = "cohort.csv"
        document = generator.make_json_dict(*arguments)
        self.assertIsNone(document["rnaseq_pipeline.fastq_index"])
        self.assertEqual("cohort", document["rnaseq_pipeline.output_report_name"])
        with self.assertRaisesRegex(ValueError, "arrays must be nonempty and aligned"):
            generator.make_json_dict(
                "human", "gencode_v39", "registry.example/rnaseq", "cohort"
            )

    def test_star_disk_type_is_opt_in_and_validated(self) -> None:
        default_document = self.make_document()
        self.assertNotIn("rnaseq_pipeline.star_disk_type", default_document)

        ssd_document = self.make_document(star_disk_type="SSD")
        self.assertEqual("SSD", ssd_document["rnaseq_pipeline.star_disk_type"])

        with self.assertRaisesRegex(ValueError, "must be HDD or SSD"):
            self.make_document(star_disk_type="LOCAL")

    def test_umi_molecule_expression_policy_requires_i1(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires a matched I1"):
            self.make_document(use_umi_molecule_expression=True)

        document = self.make_document(
            include_index=True, use_umi_molecule_expression=True
        )
        self.assertTrue(document["rnaseq_pipeline.use_umi_molecule_expression"])

        legacy_document = self.make_document(
            include_index=True, use_umi_molecule_expression=False
        )
        self.assertFalse(legacy_document["rnaseq_pipeline.use_umi_molecule_expression"])
        self.assertFalse(self.make_document()["rnaseq_pipeline.use_umi_molecule_expression"])

    def test_qc_toggles_are_default_on_and_only_disabled_values_are_emitted(self) -> None:
        default_document = self.make_document()
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
        self.assertNotIn("rnaseq_pipeline.run_multiqc", default_document)

        selective_document = self.make_document(
            run_pretrim_fastqc=False,
            run_contamination_qc=False,
            run_umi_qc=False,
        )
        for name in ("run_pretrim_fastqc", "run_contamination_qc", "run_umi_qc"):
            self.assertFalse(selective_document["rnaseq_pipeline." + name])
        for name in ("run_posttrim_fastqc", "run_alignment_qc"):
            self.assertNotIn("rnaseq_pipeline." + name, selective_document)

        multiqc_document = self.make_document(run_multiqc=True)
        self.assertTrue(multiqc_document["rnaseq_pipeline.run_multiqc"])
        for disabled_group in (
            "run_pretrim_fastqc",
            "run_posttrim_fastqc",
            "run_alignment_qc",
        ):
            with self.subTest(disabled_group=disabled_group):
                with self.assertRaisesRegex(ValueError, "MultiQC requires"):
                    self.make_document(run_multiqc=True, **{disabled_group: False})

    def test_contamination_qc_fusion_and_sampling_are_opt_in(self) -> None:
        default_document = self.make_document()
        for key in ("combine_contamination_qc", "contamination_qc_pairs"):
            self.assertNotIn("rnaseq_pipeline." + key, default_document)

        combined_document = self.make_document(combine_contamination_qc=True)
        self.assertTrue(combined_document["rnaseq_pipeline.combine_contamination_qc"])
        self.assertNotIn("rnaseq_pipeline.contamination_qc_pairs", combined_document)

        sampled_document = self.make_document(contamination_qc_pairs=1_000_000)
        self.assertEqual(1_000_000, sampled_document["rnaseq_pipeline.contamination_qc_pairs"])
        self.assertNotIn("rnaseq_pipeline.combine_contamination_qc", sampled_document)

        with self.assertRaisesRegex(ValueError, "nonnegative integer"):
            self.make_document(contamination_qc_pairs=-1)
        with self.assertRaisesRegex(ValueError, "requires run_contamination_qc"):
            self.make_document(
                run_contamination_qc=False,
                combine_contamination_qc=True,
            )

    def test_checked_runtime_profiles_match_workflow_and_only_change_resources(self) -> None:
        default_document = self.make_document()

        for filename in (
            "runtime-human-v47-small-v1.json",
            "runtime-human-v47-full-lean-v1.json",
            "runtime-human-v47-high-candidate-v1.json",
        ):
            overrides = generator.load_runtime_profile(
                REPO_ROOT / "config" / "backends" / "gcp" / filename
            )
            self.assertEqual(generator.RUNTIME_RESOURCE_KEYS, set(overrides))
            profiled_document = self.make_document(runtime_overrides=overrides)
            self.assertEqual(
                overrides,
                {
                    key: profiled_document[key]
                    for key in generator.RUNTIME_RESOURCE_KEYS
                },
            )
            self.assertEqual(
                {
                    key: value
                    for key, value in default_document.items()
                    if key not in generator.RUNTIME_RESOURCE_KEYS
                },
                {
                    key: value
                    for key, value in profiled_document.items()
                    if key not in generator.RUNTIME_RESOURCE_KEYS
                },
            )

        self.assertEqual(default_document, self.make_document())
        profile_dir = REPO_ROOT / "config" / "backends" / "gcp"
        lean = generator.load_runtime_profile(
            profile_dir / "runtime-human-v47-full-lean-v1.json"
        )
        high = generator.load_runtime_profile(
            profile_dir / "runtime-human-v47-high-candidate-v1.json"
        )

        self.assertEqual(set(lean), set(high))
        changed = {key for key in lean if lean[key] != high[key]}
        self.assertEqual({"rnaseq_pipeline.star_disk"}, changed)
        self.assertGreater(high["rnaseq_pipeline.star_disk"], lean["rnaseq_pipeline.star_disk"])
        for cpu_key in (key for key in lean if key.endswith("_ncpu")):
            ram_key = cpu_key.removesuffix("_ncpu") + "_ramGB"
            self.assertLessEqual(lean[ram_key], 8 * lean[cpu_key])
        workflow = (REPO_ROOT / "wdl/rnaseq_pipeline_scatter.wdl").read_text()
        declared = {
            "rnaseq_pipeline." + name
            for name in re.findall(
                r"\bInt\s+([A-Za-z0-9_]+_(?:ncpu|ramGB|disk))\b", workflow
            )
        }
        self.assertEqual(generator.RUNTIME_RESOURCE_KEYS, declared)

    def test_runtime_profile_rejects_incomplete_unknown_and_invalid_values(self) -> None:
        valid = {
            "schema_version": 1,
            "profile_id": "test",
            "description": "test profile",
            "overrides": {key: 1 for key in generator.RUNTIME_RESOURCE_KEYS},
        }
        cases = []

        missing = json.loads(json.dumps(valid))
        missing["overrides"].pop(next(iter(generator.RUNTIME_RESOURCE_KEYS)))
        cases.append((missing, "missing"))

        unknown = json.loads(json.dumps(valid))
        unknown["overrides"]["rnaseq_pipeline.unknown_ncpu"] = 1
        cases.append((unknown, "unknown"))

        disk_type = json.loads(json.dumps(valid))
        disk_type["overrides"]["rnaseq_pipeline.star_disk_type"] = "SSD"
        cases.append((disk_type, "unknown"))

        for invalid_value in (0, True, 1.5, "1"):
            invalid = json.loads(json.dumps(valid))
            invalid["overrides"]["rnaseq_pipeline.star_ramGB"] = invalid_value
            cases.append((invalid, "positive integers"))

        for index, (profile, message) in enumerate(cases):
            path = self.temp / "runtime-{}.json".format(index)
            path.write_text(json.dumps(profile), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, message):
                generator.load_runtime_profile(path)

    def test_invalid_profiles_and_existing_outputs_fail_before_gcs_access(self) -> None:
        path = self.temp / "invalid-runtime.json"
        path.write_text("{}", encoding="utf-8")
        arguments = self.arguments()
        arguments.runtime_profile = str(path)
        with self.assertRaisesRegex(ValueError, "schema_version"):
            generator.main(arguments)

        (self.temp / "set2_rnaseq.json").write_text("{}", encoding="utf-8")
        filesystem = FakeGcsFileSystem([], set())
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(filesystem)}):
            with self.assertRaisesRegex(ValueError, "already contains"):
                generator.main(self.arguments())
        self.assertFalse(hasattr(filesystem, "pattern"))

    def test_main_emits_explicit_runtime_profile_and_star_disk_type(self) -> None:
        r1 = "bucket/sample_R1.fastq.gz"
        r2 = "gs://bucket/sample_R2.fastq.gz"
        i1 = "gs://bucket/sample_I1.fastq.gz"
        filesystem = FakeGcsFileSystem([r1], {r2, i1})
        arguments = self.arguments()
        arguments.runtime_profile = str(
            REPO_ROOT / "config/backends/gcp/runtime-human-v47-small-v1.json"
        )
        arguments.star_disk_type = "SSD"
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(filesystem)}):
            self.assertEqual(0, generator.main(arguments))
        document = json.loads((self.temp / "set1_rnaseq.json").read_text())
        self.assertEqual(64, document["rnaseq_pipeline.star_ramGB"])
        self.assertEqual(16, document["rnaseq_pipeline.rsem_ramGB"])
        self.assertEqual(40, document["rnaseq_pipeline.markdup_ramGB"])
        self.assertEqual("SSD", document["rnaseq_pipeline.star_disk_type"])
        self.assertTrue(document["rnaseq_pipeline.use_umi_molecule_expression"])

    def test_main_checks_mates_before_writing(self) -> None:
        r1 = "bucket/sample_R1.fastq.gz"
        r2 = "gs://bucket/sample_R2.fastq.gz"
        i1 = "gs://bucket/sample_I1.fastq.gz"
        filesystem = FakeGcsFileSystem([r1], {r2, i1})
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(filesystem)}):
            self.assertEqual(0, generator.main(self.arguments()))
        document = json.loads((self.temp / "set1_rnaseq.json").read_text())
        self.assertEqual(["sample"], document["rnaseq_pipeline.sample_prefix"])
        self.assertEqual([i1], document["rnaseq_pipeline.fastq_index"])
        self.assertTrue(document["rnaseq_pipeline.use_umi_molecule_expression"])
        self.assertEqual("local-test-only", filesystem.project)

        (self.temp / "set1_rnaseq.json").unlink()
        missing = FakeGcsFileSystem([r1], set())
        with mock.patch.dict(sys.modules, {"gcsfs": self.fake_gcsfs(missing)}):
            with self.assertRaisesRegex(ValueError, "required mate/index objects are missing"):
                generator.main(self.arguments())
        self.assertFalse((self.temp / "set1_rnaseq.json").exists())

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
