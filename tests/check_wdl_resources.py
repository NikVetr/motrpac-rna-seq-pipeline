"""Evaluate actual WDL resources and execute rendered merge commands locally."""
import importlib.util
import subprocess
import tempfile
from pathlib import Path
import WDL

repo = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fixtures", repo / "tests/test_merge_results.py")
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


class Files(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        return filename


fixture = fixtures.MergeResultsTests()
fixture.setUp()
try:
    fixture.test_shuffled_files_preserve_declared_samples_dotted_ids_and_gene_order()
    for name in ("merge_results", "merge_expression"):
        task = WDL.load(str(repo / "wdl/merge_results" / (name + ".wdl"))).tasks[0]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stdlib = Files("1.0", write_dir=directory)
            inputs = {"sample_prefix": fixture.samples, "rsem_files": list(map(str, fixture.rsem.iterdir())),
                      "feature_counts_files": list(map(str, fixture.fc.iterdir())),
                      "memory": 4, "disk_space": 10, "ncpu": 1, "preemptible": 0, "docker": "unused"}
            if name == "merge_results":
                inputs.update(qc_report_files=list(map(str, fixture.qc.iterdir())), output_report_name="cohort")
            else:
                inputs["output_prefix"] = "secondary"
            env = WDL.values_from_json(inputs, task.available_inputs)
            for decl in task.postinputs:
                env = env.bind(decl.name, decl.expr.eval(env, stdlib))
            command = task.command.eval(env, stdlib).value.replace("/usr/local/src/", str(repo / "wdl/merge_results") + "/")
            subprocess.run(["bash", "-c", command], cwd=root, check=True, capture_output=True)
            for filename in ("rsem_genes_count.txt", "rsem_genes_tpm.txt", "rsem_genes_fpkm.txt", "featureCounts.txt"):
                actual = root / (("secondary_" if name == "merge_expression" else "") + filename)
                assert actual.read_bytes() == (fixture.root / filename).read_bytes()
            if name == "merge_results":
                assert (root / "cohort.csv").read_bytes() == (fixture.root / "cohort.csv").read_bytes()
            assert all(path.is_symlink() for path in (root / "rsem_files").iterdir())
            for count, floor, expected in ((1, 4, 4), (75, 4, 4), (76, 4, 8), (297, 4, 16), (600, 4, 32), (297, 64, 64)):
                env = env.bind("sample_prefix", WDL.Value.Array(WDL.Type.String(), [WDL.Value.String(str(i)) for i in range(count)]))
                env = env.bind("memory", WDL.Value.Int(floor))
                for decl in task.postinputs:
                    if decl.name != "sample_order":
                        env = env.bind(decl.name, decl.expr.eval(env, stdlib))
                assert task.runtime["memory"].eval(env, stdlib).value == f"{expected}GB"
            print(name, "rendered command, byte identity, symlinks, RAM boundaries PASS")
finally:
    fixture.tearDown()


def descendants(node):
    yield node
    for child in node.children:
        yield from descendants(child)


workflow = WDL.load(str(repo / "wdl/rnaseq_pipeline_scatter.wdl")).workflow
declarations = [node for node in descendants(workflow) if isinstance(node, WDL.Tree.Decl)
                and node.name.startswith(("inferred_molecule_rsem_", "effective_molecule_rsem_"))]
assert len(declarations) == 4
for gib, memory_floor, disk_floor, expected in (
    (0, 40, 60, (40, 60)), (12, 40, 60, (40, 60)),
    (15.21, 40, 60, (48, 71)), (30, 40, 60, (76, 130)),
    (15.21, 96, 200, (96, 200)),
):
    env = WDL.Env.Bindings().bind("molecule_rsem_input_gib", WDL.Value.Float(gib))
    env = env.bind("rsem_ramGB", WDL.Value.Int(memory_floor)).bind("rsem_disk", WDL.Value.Int(disk_floor))
    for decl in declarations:
        env = env.bind(decl.name, decl.expr.eval(env, Files("1.0")))
    assert (env.resolve("effective_molecule_rsem_memory").value,
            env.resolve("effective_molecule_rsem_scratch_gb").value) == expected
print("RSEM growth and configured floors PASS")
