from __future__ import annotations

import json
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import make_json_rnaseq as generator  # noqa: E402


class ReleaseProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def complete_profile() -> dict:
        digest = "a" * 64
        return {
            "schema_version": 1,
            "profile_id": "human_gencode_v47_test_v1",
            "organism": "human",
            "version": "gencode_v47",
            "publication_state": "published",
            "compatibility": {
                "star_index_builder": "2.7.11b",
                "star_runtime": "2.7.11b",
                "rsem_reference_builder": "1.3.3",
                "rsem_runtime": "1.3.3",
            },
            "references": {
                role: "gs://example/references/{}".format(role)
                for role in generator.REFERENCE_ROLES
            },
            "images": {
                role: "registry.example/rnaseq/{}@sha256:{}".format(role, digest)
                for role in generator.IMAGE_ROLES
            },
        }

    def write_profile(self, profile: dict, name: str = "release.json") -> Path:
        path = self.temp / name
        path.write_text(json.dumps(profile), encoding="utf-8")
        return path

    def test_legacy_v39_inputs_are_unchanged_without_a_manifest(self) -> None:
        self.assertIsNone(generator.resolve_release_inputs("human", "gencode_v39"))
        document = generator.make_json_dict(
            "human",
            "gencode_v39",
            "registry.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
        )
        expected_references = {
            "rnaseq_pipeline.star_index": "gs://omicspipelines-public-resources/rnaseq/references/human/hg38_v39_star_index.tar.gz",
            "rnaseq_pipeline.gtf_file": "gs://omicspipelines-public-resources/rnaseq/references/human/GRCh38.v39.primary_assembly.annotation.gtf",
            "rnaseq_pipeline.rsem_reference": "gs://omicspipelines-public-resources/rnaseq/references/human/hg38_rsem_reference.tar.gz",
            "rnaseq_pipeline.globin_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/hs_globin.tar.gz",
            "rnaseq_pipeline.rrna_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/hs_rRNA.tar.gz",
            "rnaseq_pipeline.phix_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/phix.tar.gz",
            "rnaseq_pipeline.ref_flat": "gs://omicspipelines-public-resources/rnaseq/references/human/refFlat_hg38_v39.txt",
        }
        legacy_image_names = {
            "fastqc_docker": "fastqc",
            "attach_umi_docker": "umi_attach",
            "cutadapt_docker": "cutadapt",
            "star_docker": "star",
            "feature_counts_docker": "feature_counts",
            "rsem_docker": "rsem",
            "bowtie_docker": "bowtie",
            "picard_docker": "picard",
            "umi_dup_docker": "umi_dup",
            "samtools_docker": "samtools",
            "collect_qc_docker": "collect_qc",
            "merge_results_docker": "merge_results",
        }
        expected_images = {
            "rnaseq_pipeline.{}".format(role): "registry.example/rnaseq/{}:latest".format(
                image_name
            )
            for role, image_name in legacy_image_names.items()
        }
        self.assertEqual(
            expected_references,
            {key: document[key] for key in expected_references},
        )
        self.assertEqual(
            expected_images,
            {key: document[key] for key in expected_images},
        )

    def test_complete_manifest_maps_all_reference_and_image_inputs(self) -> None:
        profile = self.complete_profile()
        inputs = generator.load_release_manifest(
            self.write_profile(profile), "human", "gencode_v47"
        )
        expected_roles = generator.REFERENCE_ROLES | generator.IMAGE_ROLES
        self.assertEqual(
            {"rnaseq_pipeline.{}".format(role) for role in expected_roles},
            set(inputs),
        )

        document = generator.make_json_dict(
            "human",
            "gencode_v47",
            "ignored.example/rnaseq",
            "cohort",
            ["gs://example/sample_R1.fastq.gz"],
            ["gs://example/sample_R2.fastq.gz"],
            None,
            ["sample"],
            inputs,
        )
        for section_name in ("references", "images"):
            for role, value in profile[section_name].items():
                self.assertEqual(value, document["rnaseq_pipeline." + role])

    def test_published_builtin_v47_is_complete_and_immutable(self) -> None:
        self.assertIn(("human", "gencode_v47"), generator.SUPPORTED_REFERENCES)
        inputs = generator.resolve_release_inputs("human", "gencode_v47")
        self.assertEqual(
            {
                "rnaseq_pipeline.{}".format(role)
                for role in generator.REFERENCE_ROLES | generator.IMAGE_ROLES
            },
            set(inputs),
        )
        for role in generator.IMAGE_ROLES:
            self.assertRegex(inputs["rnaseq_pipeline." + role], r"@sha256:[0-9a-f]{64}$")
        for role in generator.REFERENCE_ROLES:
            self.assertTrue(inputs["rnaseq_pipeline." + role].startswith("gs://"))

    def test_manifest_identity_and_compatibility_are_enforced(self) -> None:
        cases = []

        wrong_schema = self.complete_profile()
        wrong_schema["schema_version"] = 2
        cases.append((wrong_schema, "schema_version"))

        wrong_identity = self.complete_profile()
        wrong_identity["organism"] = "rat"
        cases.append((wrong_identity, "declares rat"))

        missing_reference = self.complete_profile()
        del missing_reference["references"]["star_index"]
        cases.append((missing_reference, "missing star_index"))

        empty_reference = self.complete_profile()
        empty_reference["references"]["star_index"] = ""
        cases.append((empty_reference, "empty execution values"))

        star_mismatch = self.complete_profile()
        star_mismatch["compatibility"]["star_index_builder"] = "2.7.0d"
        cases.append((star_mismatch, "STAR index-builder/runtime"))

        rsem_mismatch = self.complete_profile()
        rsem_mismatch["compatibility"]["rsem_reference_builder"] = "1.3.0"
        cases.append((rsem_mismatch, "RSEM reference-builder/runtime"))

        mutable_image = self.complete_profile()
        mutable_image["images"]["star_docker"] = "registry.example/star:latest"
        cases.append((mutable_image, "immutable sha256"))

        for index, (profile, message) in enumerate(cases):
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    generator.load_release_manifest(
                        self.write_profile(deepcopy(profile), "case{}.json".format(index)),
                        "human",
                        "gencode_v47",
                    )


if __name__ == "__main__":
    unittest.main()
