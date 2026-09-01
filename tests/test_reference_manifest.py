import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = REPO_ROOT / "config/references/human-grch38-gencode-v47.json"


class ReferenceManifestTests(unittest.TestCase):
    def test_v47_bundle_is_annotation_baked_and_checksum_locked(self):
        manifest = json.loads(MANIFEST.read_text())
        self.assertEqual("GENCODE v47 (Ensembl 113)", manifest["annotation"]["release"])
        self.assertEqual(
            "comprehensive primary-assembly annotation",
            manifest["annotation"]["scope"],
        )
        self.assertEqual(78932, manifest["annotation"]["genes"])
        self.assertEqual(387944, manifest["annotation"]["transcripts"])
        self.assertTrue(manifest["artifacts"]["star"]["builder"]["annotation_baked"])
        self.assertEqual(100, manifest["artifacts"]["star"]["builder"]["sjdbOverhang"])

        expected_artifacts = {
            "star": (
                "2f0129bce1583341ea0ed97c700c22740c98f2681fb0a1d0f63d88188901c6a5",
                27102859416,
                "2.7.11b",
            ),
            "rsem": (
                "db190d069ed575ca48379931ae8f605ec5d87377e9a87950d30e38d57c669ddc",
                312996761,
                "1.3.3",
            ),
            "refFlat": (
                "1d985c67b63b8ae93118222ba89277634d7491230294a344edcddc318d447459",
                72603520,
                "kent-v479",
            ),
        }
        for name, (sha256, size, version) in expected_artifacts.items():
            artifact = manifest["artifacts"][name]
            self.assertEqual(sha256, artifact["sha256"])
            self.assertEqual(size, artifact["bytes"])
            self.assertEqual(version, artifact["builder"]["version"])
            uri = artifact["distribution_uri"]
            self.assertTrue(uri.startswith("https://"))
            self.assertIn("sha256-{}".format(sha256), uri)
            self.assertTrue(artifact["gcs_uri"].startswith("gs://"))
            self.assertRegex(artifact["gcs_generation"], r"^[1-9][0-9]+$")

        annotation = manifest["annotation"]
        self.assertIn(
            "sha256-{}".format(annotation["decompressed_sha256"]),
            annotation["distribution_uri"],
        )
        self.assertRegex(annotation["gcs_generation"], r"^[1-9][0-9]+$")
        self.assertIn("published", manifest["publication_state"])

        expected_sha256 = {
            "annotation": (
                "f02ee3e1c8e7fd9be264be6d0b974feb225a1e9d6c81915ee271804b060bd8c0",
                "7478f1c14e4915c32822e03cc7255a52cf0ae5ff4d08e926774dbb96b48b0c9b",
            ),
            "genome": (
                "fca1b272425c11d5fa3ab11b2052e815d1b00de722b49617d7819cc998eb1fc1",
                "e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b",
            ),
        }
        for name, (compressed_sha256, decompressed_sha256) in expected_sha256.items():
            source = manifest[name]
            self.assertRegex(source["source_uri"], r"^https://")
            self.assertRegex(source["compressed_md5"], r"^[0-9a-f]{32}$")
            self.assertEqual(compressed_sha256, source["compressed_sha256"])
            self.assertEqual(decompressed_sha256, source["decompressed_sha256"])


if __name__ == "__main__":
    unittest.main()
