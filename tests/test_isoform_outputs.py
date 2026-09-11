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
            inputs = root / "inputs"
            inputs.mkdir()
            order = root / "order"
            order.write_text("b\na\n")
            header = "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n"
            a = "tx.1\tg.1\t100\t70.5\t1.20\t2.30\t3.40\t50\ntx.2\tg.1\t200\t170.5\t0.00\t0.00\t0.00\t50\n"
            b = a.replace("1.20", "4.50")
            (inputs / "a.isoforms.results").write_text(header + a)
            target = inputs / "b.isoforms.results"
            target.write_text(header + b)
            merger.merge_isoforms(inputs, order, root)
            self.assertEqual("transcript_id\tb\ta\ntx.1\t4.50\t1.20\ntx.2\t0.00\t0.00\n",
                             (root / "rsem_isoforms_count.txt").read_text())
            self.assertEqual(header + a, (inputs / "a.isoforms.results").read_text())
            for index, broken in enumerate((b.replace("g.1", "g.2"), b.splitlines()[0] + "\n",
                                             b + b.splitlines()[0] + "\n", b.replace("tx.2", "tx.1"))):
                output = root / str(index)
                output.mkdir()
                target.write_text(header + broken)
                with self.assertRaises(ValueError):
                    merger.merge_isoforms(inputs, order, output)

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
