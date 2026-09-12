"""Evaluate the production WDL's per-sample policy and metadata on mixed inputs."""
import csv
import argparse
import json
import subprocess
import tempfile
from pathlib import Path
import WDL

ROOT = Path(__file__).resolve().parents[1]
workflow = WDL.load(str(ROOT / "wdl/rnaseq_pipeline_scatter.wdl")).workflow
merge = next(node for node in workflow.body if isinstance(node, WDL.Tree.Call) and node.name == "merge_results")
metadata_writer = next(part.expr for part in merge.callee.command.parts
                       if isinstance(part, WDL.Expr.Placeholder) and str(part.expr) == "write_tsv(expression_metadata_rows)")


def descendants(node):
    yield node
    for child in node.children:
        yield from descendants(child)


declarations = {node.name: node for node in descendants(workflow) if isinstance(node, WDL.Tree.Decl)}
global_names = ("index_uris", "has_fastq_index", "index_length_contract", "index_length_valid",
                "expression_policy_contract", "expression_policy_valid")
sample_names = ("index_uri", "sample_index", "umi_expression_input_contract", "umi_expression_inputs_valid",
                "use_index_reads", "use_sample_umi_expression", "run_all_read_expression",
                "expression_mode", "umi_status", "expression_metadata_row")
with tempfile.TemporaryDirectory() as directory:
    stdlib = WDL.StdLib.Base("1.0", write_dir=directory)
    for indexes, allow, molecule, retain, expected in (
        (["a.I1", "b.I1"], False, True, False, ["umi_molecules", "umi_molecules"]),
        (["a.I1", ""], True, True, False, ["umi_molecules", "all_read"]),
        (["", ""], True, True, False, ["all_read", "all_read"]),
        (None, True, True, False, ["all_read", "all_read"]),
        (None, False, False, False, ["all_read", "all_read"]),
        (["a.I1", ""], True, True, True, ["umi_molecules", "all_read"]),
        (["a.I1", ""], False, True, False, None),
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
            rows = []
            for index in range(2):
                sample_env = env.bind("i", WDL.Value.Int(index))
                for name in sample_names:
                    sample_env = sample_env.bind(name, declarations[name].expr.eval(sample_env, stdlib))
                mode = sample_env.resolve("expression_mode").value
                if expected is not None:
                    assert mode == expected[index]
                assert sample_env.resolve("run_all_read_expression").value == (retain or mode == "all_read")
                row = sample_env.resolve("expression_metadata_row")
                assert row.json[3] == ("0" if mode == "umi_molecules" else "1")
                rows.append(row)
        except WDL.Error.EvalError:
            assert expected is None
            continue
        assert expected is not None
        metadata_env = env.bind("expression_metadata_row", WDL.Value.Array(WDL.Type.Array(WDL.Type.String()), rows))
        metadata_env = metadata_env.bind("expression_metadata_rows",
            merge.inputs["expression_metadata_rows"].eval(metadata_env, stdlib))
        # Host paths are used only for this local expression check.
        stdlib._virtualize_filename = lambda filename: filename
        path = metadata_writer.eval(metadata_env, stdlib).value
        with open(path) as handle:
            actual = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["sample"] for row in actual] == ["a", "b"]
        assert [row["expression_mode"] for row in actual] == expected
print("Mixed, absent, strict, all-read and retained-secondary WDL policies PASS")

# Static validation does not expose Cromwell's nullable-array evaluation bugs.
parser = argparse.ArgumentParser()
parser.add_argument("--cromwell-jar", type=Path)
parser.add_argument("--java", default="java")
args = parser.parse_args()
if args.cromwell_jar:
    def declaration(name):
        decl = declarations[name]
        return f"{decl.type} {name}" + (f" = {decl.expr}" if decl.expr else "")

    probe = "version 1.0\nworkflow policy {\ninput {\n" + "\n".join(declaration(name) for name in inputs) + "\n}\n"
    probe += "\n".join(declaration(name) for name in global_names)
    probe += "\nscatter (i in range(length(fastq1))) {\n" + "\n".join(declaration(name) for name in sample_names)
    probe += "\n}\noutput { Array[String] modes = expression_mode }\n}\n"
    for indexes, allow, expected in ((["/etc/hosts", ""], True, ["umi_molecules", "all_read"]),
                                     (None, True, ["all_read", "all_read"]),
                                     (["/etc/hosts", ""], False, None)):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "policy.wdl").write_text(probe)
            values = dict(inputs, fastq1=["/etc/hosts"] * 2, fastq2=["/etc/hosts"] * 2,
                          fastq_index=indexes, allow_missing_umis=allow)
            (root / "input.json").write_text(json.dumps({"policy." + key: value for key, value in values.items()}))
            command = [args.java, "-Xmx1g", "-Dbackend.default=Local", "-jar", str(args.cromwell_jar.resolve()),
                       "run", "policy.wdl", "-i", "input.json", "-m", "metadata.json"]
            completed = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=120)
            metadata = json.loads((root / "metadata.json").read_text())
            if expected is None:
                assert metadata["status"] == "Failed" and "umi_expression_inputs_valid" in json.dumps(metadata["failures"])
            else:
                assert completed.returncode == 0, completed.stdout + completed.stderr
                assert metadata["outputs"]["policy.modes"] == expected
    print("Cromwell mixed, absent, and strict expression evaluation PASS (no worker tasks)")
