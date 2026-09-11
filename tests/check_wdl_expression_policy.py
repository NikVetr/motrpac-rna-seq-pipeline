"""Evaluate the production WDL's per-sample policy and metadata on mixed inputs."""
import csv
import tempfile
from pathlib import Path
import WDL

ROOT = Path(__file__).resolve().parents[1]
workflow = WDL.load(str(ROOT / "wdl/rnaseq_pipeline_scatter.wdl")).workflow


def descendants(node):
    yield node
    for child in node.children:
        yield from descendants(child)


declarations = {node.name: node for node in descendants(workflow) if isinstance(node, WDL.Tree.Decl)}
global_names = ("has_fastq_index", "index_length_contract", "index_length_valid",
                "umi_expression_input_contract", "umi_expression_inputs_valid",
                "expression_policy_contract", "expression_policy_valid")
sample_names = ("sample_index", "use_index_reads", "use_sample_umi_expression", "run_all_read_expression",
                "expression_mode", "umi_status", "expression_metadata_row")
with tempfile.TemporaryDirectory() as directory:
    stdlib = WDL.StdLib.Base("1.0", write_dir=directory)
    for indexes, allow, molecule, retain, expected in (
        (["a.I1", "b.I1"], False, True, False, ["umi_molecules", "umi_molecules"]),
        (["a.I1", None], True, True, False, ["umi_molecules", "all_read"]),
        ([None, None], True, True, False, ["all_read", "all_read"]),
        (None, True, True, False, ["all_read", "all_read"]),
        (None, False, False, False, ["all_read", "all_read"]),
        (["a.I1", None], True, True, True, ["umi_molecules", "all_read"]),
        (["a.I1", None], False, True, False, None),
        (None, False, True, False, None),
        (["a.I1"], True, True, False, None),
    ):
        inputs = {"fastq1": ["a.R1", "b.R1"], "fastq2": ["a.R2", "b.R2"],
                  "sample_prefix": ["a", "b"], "fastq_index": indexes, "allow_missing_umis": allow,
                  "use_umi_molecule_expression": molecule, "retain_all_read_expression": retain,
                  "reference_release": "gencode_v50"}
        env = WDL.values_from_json(inputs, workflow.available_inputs)
        try:
            for name in global_names:
                env = env.bind(name, declarations[name].expr.eval(env, stdlib))
        except WDL.Error.EvalError:
            assert expected is None
            continue
        assert expected is not None
        rows = []
        for index, mode in enumerate(expected):
            sample_env = env.bind("i", WDL.Value.Int(index))
            for name in sample_names:
                sample_env = sample_env.bind(name, declarations[name].expr.eval(sample_env, stdlib))
            assert sample_env.resolve("expression_mode").value == mode
            assert sample_env.resolve("run_all_read_expression").value == (retain or mode == "all_read")
            row = sample_env.resolve("expression_metadata_row")
            assert row.json[3] == ("0" if mode == "umi_molecules" else "1")
            rows.append(row)
        metadata_env = env.bind("expression_metadata_row", WDL.Value.Array(WDL.Type.Array(WDL.Type.String()), rows))
        # Host paths are used only for this local expression check.
        stdlib._virtualize_filename = lambda filename: filename
        path = declarations["expression_metadata"].expr.eval(metadata_env, stdlib).value
        with open(path) as handle:
            actual = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["sample"] for row in actual] == ["a", "b"]
        assert [row["expression_mode"] for row in actual] == expected
print("Mixed, absent, strict, all-read and retained-secondary WDL policies PASS")
