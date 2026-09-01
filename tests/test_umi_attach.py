from __future__ import annotations

import gzip
import shutil
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
ATTACH_SCRIPT = REPO_ROOT / "wdl" / "attach_umi" / "UMI_attach.awk"
AWK = shutil.which("gawk") or shutil.which("awk")


def fastq_record(name: str, sequence: str, quality: str | None = None) -> bytes:
    quality = quality if quality is not None else "I" * len(sequence)
    return f"@{name} read metadata\n{sequence}\n+\n{quality}\n".encode()


class AttachUmiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def run_attach(self, biological: bytes, index: bytes) -> subprocess.CompletedProcess:
        with gzip.open(self.temp / "index.fastq.gz", "wb") as handle:
            handle.write(index)
        return subprocess.run(
            [AWK, "-v", "Ifq=index.fastq.gz", "-f", str(ATTACH_SCRIPT)],
            cwd=self.temp,
            input=biological,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def run_task_shape(
        self,
        biological_r1: bytes,
        biological_r2: bytes,
        index: bytes,
        corrupt_r2_gzip: bool = False,
    ) -> subprocess.CompletedProcess:
        r1 = self.temp / "input_R1.fastq.gz"
        r2 = self.temp / "input_R2.fastq.gz"
        i1 = self.temp / "input_I1.fastq.gz"
        for path, content in ((r1, biological_r1), (r2, biological_r2), (i1, index)):
            with gzip.open(path, "wb") as handle:
                handle.write(content)
        if corrupt_r2_gzip:
            r2.write_bytes(b"not a gzip stream")

        quoted = {name: shlex.quote(str(value)) for name, value in {
            "awk": AWK,
            "script": ATTACH_SCRIPT,
            "r1": r1,
            "r2": r2,
            "i1": i1,
        }.items()}
        command = """
set -euo pipefail
mkdir fastq_attach
r1_tmp=fastq_attach/sample_R1.fastq.gz.tmp
r2_tmp=fastq_attach/sample_R2.fastq.gz.tmp
r1_pid=
r2_pid=
trap 'rm -f -- "$r1_tmp" "$r2_tmp"; kill $r1_pid $r2_pid 2>/dev/null || true' EXIT
(
set -euo pipefail
gzip -cd -- {r1} | {awk} -v Ifq={i1} -f {script} | gzip -c > "$r1_tmp"
) &
r1_pid=$!
(
set -euo pipefail
gzip -cd -- {r2} | {awk} -v Ifq={i1} -f {script} | gzip -c > "$r2_tmp"
) &
r2_pid=$!
set +e
wait "$r1_pid"; r1_status=$?
wait "$r2_pid"; r2_status=$?
set -e
(( r1_status == 0 && r2_status == 0 ))
gzip -t "$r1_tmp"
gzip -t "$r2_tmp"
mv -- "$r1_tmp" fastq_attach/sample_R1.fastq.gz
mv -- "$r2_tmp" fastq_attach/sample_R2.fastq.gz
trap - EXIT
""".format(**quoted)
        return subprocess.run(
            ["bash", "-c", command],
            cwd=self.temp,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_valid_records_preserve_legacy_header_format(self) -> None:
        result = self.run_attach(
            fastq_record("read1", "ACGT"),
            fastq_record("read1", "AACCGGTT"),
        )
        self.assertEqual(0, result.returncode, result.stderr.decode())
        self.assertEqual(
            b"@read1:AACCGGTT read metadata\nACGT\n+\nIIII\n",
            result.stdout,
        )

        no_metadata = self.run_attach(
            b"@read1\nACGT\n+\nIIII\n",
            b"@read1\nAACCGGTT\n+\nIIIIIIII\n",
        )
        self.assertEqual(0, no_metadata.returncode, no_metadata.stderr.decode())
        self.assertEqual(b"@read1:AACCGGTT\nACGT\n+\nIIII\n", no_metadata.stdout)

    def test_rejects_missing_extra_or_mismatched_index_records(self) -> None:
        two_biological = fastq_record("read1", "ACGT") + fastq_record("read2", "TGCA")
        one_index = fastq_record("read1", "AACCGGTT")
        self.assertNotEqual(0, self.run_attach(two_biological, one_index).returncode)

        one_biological = fastq_record("read1", "ACGT")
        two_index = one_index + fastq_record("read2", "TTGGCCAA")
        self.assertNotEqual(0, self.run_attach(one_biological, two_index).returncode)

        mismatched = fastq_record("different", "AACCGGTT")
        self.assertNotEqual(0, self.run_attach(one_biological, mismatched).returncode)

    def test_rejects_truncated_fastq_and_non_eight_base_umi(self) -> None:
        truncated = b"@read1 read metadata\nACGT\n+\n"
        valid_index = fastq_record("read1", "AACCGGTT")
        self.assertNotEqual(0, self.run_attach(truncated, valid_index).returncode)

        short_index = fastq_record("read1", "AACCGGT")
        self.assertNotEqual(
            0, self.run_attach(fastq_record("read1", "ACGT"), short_index).returncode
        )

    def test_rejects_corrupt_index_gzip(self) -> None:
        index_path = self.temp / "index.fastq.gz"
        index_path.write_bytes(b"not a gzip stream")
        result = subprocess.run(
            [AWK, "-v", "Ifq=index.fastq.gz", "-f", str(ATTACH_SCRIPT)],
            cwd=self.temp,
            input=fastq_record("read1", "ACGT"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertNotEqual(0, result.returncode)

    def test_task_shape_publishes_both_mates_only_after_validation(self) -> None:
        biological = fastq_record("read1", "ACGT")
        index = fastq_record("read1", "AACCGGTT")
        valid = self.run_task_shape(biological, biological, index)
        self.assertEqual(0, valid.returncode, valid.stderr.decode())
        self.assertTrue((self.temp / "fastq_attach/sample_R1.fastq.gz").is_file())
        self.assertTrue((self.temp / "fastq_attach/sample_R2.fastq.gz").is_file())

        shutil.rmtree(self.temp / "fastq_attach")
        mismatched_r2 = fastq_record("different", "ACGT")
        failed = self.run_task_shape(biological, mismatched_r2, index)
        self.assertNotEqual(0, failed.returncode)
        self.assertFalse((self.temp / "fastq_attach/sample_R1.fastq.gz").exists())
        self.assertFalse((self.temp / "fastq_attach/sample_R2.fastq.gz").exists())

        shutil.rmtree(self.temp / "fastq_attach")
        corrupt = self.run_task_shape(
            biological, biological, index, corrupt_r2_gzip=True
        )
        self.assertNotEqual(0, corrupt.returncode)
        self.assertFalse((self.temp / "fastq_attach/sample_R1.fastq.gz").exists())
        self.assertFalse((self.temp / "fastq_attach/sample_R2.fastq.gz").exists())

    def test_wdl_and_image_use_the_validated_repository_script(self) -> None:
        wdl = (REPO_ROOT / "wdl/attach_umi/attach_umi.wdl").read_text()
        dockerfile = (REPO_ROOT / "dockerfiles/umi_attach.Dockerfile").read_text()
        self.assertIn("set -euo pipefail", wdl)
        self.assertIn('wait "$r1_pid"', wdl)
        self.assertIn('wait "$r2_pid"', wdl)
        self.assertIn("UMI attachment mate failures", wdl)
        self.assertLess(wdl.index('gzip -t "$r2_tmp"'), wdl.index('mv -- "$r1_tmp"'))
        self.assertIn(
            "COPY wdl/attach_umi/UMI_attach.awk /usr/local/src/UMI_attach.awk",
            dockerfile,
        )


if __name__ == "__main__":
    unittest.main()
