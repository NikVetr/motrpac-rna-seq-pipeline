import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class RatReferenceTests(unittest.TestCase):
    def test_contig_conversion_preserves_scaffolds_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta, gtf = root / "genome.fa", root / "genes.gtf"
            names = ["1", "20", "X", "Y", "MT", "MU150001.1", "chr2"]
            expected = ["chr1", "chr20", "chrX", "chrY", "chrM", "MU150001.1", "chr2"]
            fasta.write_text("".join(">{}{}\nACGT\n".format(name, " description" if i % 2 else "")
                                     for i, name in enumerate(names)))
            gtf.write_text("# annotation\n" + "".join(
                '{}\tEnsembl\texon\t1\t4\t.\t+\t.\tgene_id "g{}"; transcript_id "t{}";\n'.format(name, i, i)
                for i, name in enumerate(names)))
            command = ["bash", str(ROOT / "scripts/add_chr_prefix.sh"), str(fasta), str(gtf)]
            for _ in range(2):
                subprocess.run(command, check=True, capture_output=True, text=True)
                self.assertEqual(expected, [line[1:].split()[0] for line in fasta.read_text().splitlines() if line.startswith(">")])
                self.assertEqual(expected, [line.split("\t")[0] for line in gtf.read_text().splitlines() if not line.startswith("#")])

    def test_mapped_categories_agree_for_ensembl_and_ucsc_names(self):
        wdl = (ROOT / "wdl/compute_mapped/mapped.wdl").read_text()
        program = wdl.split("awk -v name=~{SID} '\n", 1)[1].split("\n        ' ~{SID}_aligned_chr_info.txt", 1)[0]
        for names in (["1", "20", "X", "Y", "MT", "MU150001.1"],
                      ["chr1", "chr20", "chrX", "chrY", "chrM", "MU150001.1"]):
            records = "".join("{}\t100\t{}\t0\n".format(name, count)
                              for name, count in zip(names, [30, 20, 10, 5, 15, 20]))
            result = subprocess.run(["awk", "-v", "name=rat", program], input=records,
                                    text=True, capture_output=True, check=True)
            self.assertEqual(["rat", "10", "5", "15", "50", "20"], result.stdout.splitlines()[1].split("\t"))


if __name__ == "__main__":
    unittest.main()
