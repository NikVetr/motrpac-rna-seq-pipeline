import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class ModernToolTests(unittest.TestCase):
    def test_v47_tool_profile_is_immutable_and_matches_dockerfiles(self):
        manifest = json.loads(
            (REPO_ROOT / "config/tool-images-gencode-v47.json").read_text()
        )
        self.assertEqual("human_gencode_v47", manifest["profile"])
        self.assertEqual(
            {
                "bowtie2", "cutadapt", "fastqc", "picard", "rsem",
                "multiqc", "samtools", "star", "subread", "umi_tools",
            },
            set(manifest["images"]),
        )
        for record in manifest["images"].values():
            self.assertRegex(record["uri"], r"@sha256:[0-9a-f]{64}$")

        dockerfiles = {
            "bowtie2": "bowtie.Dockerfile",
            "cutadapt": "cutadapt.Dockerfile",
            "fastqc": "fastqc.Dockerfile",
            "picard": "picard.Dockerfile",
            "rsem": "rsem.Dockerfile",
            "samtools": "samtools.Dockerfile",
            "star": "star.Dockerfile",
            "subread": "feature_counts.Dockerfile",
        }
        for tool, filename in dockerfiles.items():
            content = (REPO_ROOT / "dockerfiles" / filename).read_text().strip()
            self.assertEqual("FROM " + manifest["images"][tool]["uri"], content)

    def test_updated_command_interfaces_are_explicit(self):
        cutadapt = (REPO_ROOT / "wdl/cutadapt/cutadapt.wdl").read_text()
        self.assertIn("--cores ~{ncpu}", cutadapt)
        self.assertIn("--pair-filter any", cutadapt)
        self.assertIn("--compression-level 1", cutadapt)

        markdup = (REPO_ROOT / "wdl/mark_duplicates/mark_duplicates.wdl").read_text()
        metrics = (
            REPO_ROOT / "wdl/collect_rnaseq_metrics/collect_rnaseq_metrics.wdl"
        ).read_text()
        self.assertIn("picard -Xmx32g MarkDuplicates", markdup)
        self.assertIn("picard -Xmx~{memory}g CollectRnaSeqMetrics", metrics)
        self.assertNotIn("picard.jar", markdup + metrics)

    def test_counting_is_forward_stranded_and_fragment_based(self):
        feature_counts = (REPO_ROOT / "wdl/feature_counts/fc.wdl").read_text()
        for argument in ("-T ~{ncpu}", "--countReadPairs", "-s 1"):
            self.assertIn(argument, feature_counts)

        rsem = (REPO_ROOT / "wdl/rsem_exp/rsem.wdl").read_text()
        self.assertIn("--paired-end", rsem)
        self.assertIn("--forward-prob 1", rsem)
        self.assertNotIn("--forward-prob 0.5", rsem)


if __name__ == "__main__":
    unittest.main()
