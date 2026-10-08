import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("merge_rsem", ROOT / "wdl/merge_results/merge_rsem.py")
merger = importlib.util.module_from_spec(spec)
spec.loader.exec_module(merger)


class IsoformOutputTests(unittest.TestCase):
    def test_streaming_merge_preserves_values_and_rejects_incompatible_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs, genes = root / "inputs", root / "genes"
            inputs.mkdir()
            genes.mkdir()
            order = root / "order"
            order.write_text("b\na\n")
            header = "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n"
            a = "tx.1\tg.1\t100\t70.5\t1.20\t1000000.00\t3.40\t100\ntx.2\tg.1\t200\t170.5\t0.00\t0.00\t0.00\t0\n"
            b = a.replace("1.20", "4.50")
            gene_header = "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n"
            for sample, count in (("a", "1.20"), ("b", "4.50")):
                (genes / f"{sample}.genes.results").write_text(
                    gene_header + f"g.1\ttx.1,tx.2\t100.00\t70.50\t{count}\t1000000.00\t3.40\n")
            (inputs / "a.isoforms.results").write_text(header + a)
            target = inputs / "b.isoforms.results"
            target.write_text(header + b)
            merger.merge_isoforms(inputs, order, root, genes)
            self.assertEqual("transcript_id\tb\ta\ntx.1\t4.50\t1.20\ntx.2\t0.00\t0.00\n",
                             (root / "rsem_isoforms_count.txt").read_text())
            self.assertEqual("transcript_id\tb\ta\ntx.1\t70.5\t70.5\ntx.2\t170.5\t170.5\n",
                             (root / "rsem_isoforms_effective_length.txt").read_text())
            self.assertEqual("transcript_id\tgene_id\tlength\ntx.1\tg.1\t100\ntx.2\tg.1\t200\n",
                             (root / "rsem_transcripts.tsv").read_text())
            self.assertEqual(header + a, (inputs / "a.isoforms.results").read_text())
            # Incompatible rows, and a gene count that differs from its isoform sum, fail loudly.
            for index, broken in enumerate((b.replace("g.1", "g.2"), b.splitlines()[0] + "\n",
                                             b + b.splitlines()[0] + "\n", b.replace("tx.2", "tx.1"),
                                             b.replace("4.50", "5.50"))):
                output = root / str(index)
                output.mkdir()
                target.write_text(header + broken)
                with self.assertRaises(ValueError):
                    merger.merge_isoforms(inputs, order, output, genes)

    def test_expression_covariate_joins_in_matrix_order(self):
        spec = importlib.util.spec_from_file_location("metadata", ROOT / "scripts/prepare_sample_metadata.py")
        metadata = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(metadata)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = {name: root / name for name in ("matrix", "study", "qc", "expression", "output")}
            paths["matrix"].write_text("transcript_id\tb\ta\ntx\t1\t2\n")
            paths["study"].write_text("sample,pid\na,001\nb,002\n")
            paths["qc"].write_text("sample,pct_umi_dup\na,10\nb,\n")
            paths["expression"].write_text("sample\tnot_deduplicated\na\t0\nb\t1\n")
            metadata.prepare(paths["matrix"], paths["study"], paths["qc"], paths["output"],
                             ["pid"], paths["expression"])
            with paths["output"].open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([("b", "1"), ("a", "0")], [(row["sample"], row["not_deduplicated"]) for row in rows])
            paths["expression"].write_text("sample\tnot_deduplicated\na\t0\n")
            with self.assertRaisesRegex(ValueError, "samples differ"):
                metadata.prepare(paths["matrix"], paths["study"], paths["qc"], root / "bad", ["pid"], paths["expression"])


if __name__ == "__main__":
    unittest.main()
