from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


def source(relative_path):
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


class WdlIoContractTests(unittest.TestCase):
    def test_bowtie_discards_only_unused_alignment_records(self):
        wdl = source("wdl/bowtie2_align/bowtie2_align.wdl")
        self.assertIn("set -euo pipefail", wdl)
        self.assertIn("--local -S /dev/null", wdl)
        self.assertIn("invalid Bowtie2 terminal alignment percentage", wdl)
        self.assertNotIn("File bowtie2_output", wdl)
        self.assertIn("File bowtie2_log", wdl)
        self.assertIn("File bowtie2_report", wdl)

    def test_markduplicates_publishes_metrics_only(self):
        wdl = source("wdl/mark_duplicates/mark_duplicates.wdl")
        self.assertIn("ln -s /dev/null ~{output_bam}", wdl)
        self.assertIn("picard -Xmx32g MarkDuplicates", wdl)
        self.assertIn("CREATE_INDEX=false", wdl)
        self.assertIn("File metrics", wdl)
        self.assertNotIn("File bam_file", wdl)
        self.assertNotIn("File bam_index", wdl)

    def test_chromosome_summary_streams_primary_records(self):
        wdl = source("wdl/compute_mapped/mapped.wdl")
        self.assertIn("samtools view -H", wdl)
        self.assertIn("samtools view -F 0x900", wdl)
        self.assertNotIn("samtools view -b", wdl)
        self.assertNotIn("samtools index", wdl)
        self.assertIn('File aligned_chrinfo = "${SID}_aligned_chr_info.txt"', wdl)

    def test_fastqc_collects_both_parallel_statuses(self):
        wdl = source("wdl/fastqc/fastqc.wdl")
        self.assertIn("set -euo pipefail", wdl)
        self.assertIn("r1_pid=$!", wdl)
        self.assertIn("r2_pid=$!", wdl)
        self.assertIn('wait "$r1_pid"', wdl)
        self.assertIn('wait "$r2_pid"', wdl)
        self.assertIn("exactly one ZIP and HTML report per mate", wdl)

    def test_star_emits_truthful_sample_read_group_without_unused_index(self):
        wdl = source("wdl/star_align/star.wdl")
        self.assertIn(
            "--outSAMattrRGline ID:~{prefix} SM:~{prefix} PL:ILLUMINA", wdl
        )
        self.assertNotIn("samtools index", wdl)
        self.assertNotIn("File bam_index", wdl)

    def test_multiqc_does_not_remove_localized_inputs(self):
        for path in (
            "wdl/multiqc/multiqc.wdl",
            "wdl/multiqc/multiqc_postalign.wdl",
        ):
            self.assertNotIn("rm $FILE", source(path))


if __name__ == "__main__":
    unittest.main()
