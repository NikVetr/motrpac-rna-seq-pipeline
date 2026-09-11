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
                inputs["expression_metadata_rows"] = [["sample", "not_deduplicated"]] + [[sample, "0"] for sample in fixture.samples]
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
                metadata = task.outputs[-1].expr.eval(env, stdlib).value
                assert (root / metadata).read_text() == "sample\tnot_deduplicated\n" + "".join(sample + "\t0\n" for sample in fixture.samples)
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

task = WDL.load(str(repo / "wdl/merge_results/merge_isoforms.wdl")).tasks[0]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "source"
    source.mkdir()
    for sample, count in (("a", "1.20"), ("b", "2.30")):
        (source / (sample + ".isoforms.results")).write_text(
            "transcript_id\tgene_id\texpected_count\tTPM\tFPKM\n"
            + "tx.1\tg.1\t" + count + "\t3.40\t5.60\n")
    stdlib = Files("1.0", write_dir=directory)
    env = WDL.values_from_json({"sample_prefix": ["b", "a"],
        "rsem_files": list(map(str, source.iterdir())), "memory": 1,
        "disk_space": 1, "ncpu": 1, "preemptible": 0, "docker": "unused"}, task.available_inputs)
    for decl in task.postinputs:
        env = env.bind(decl.name, decl.expr.eval(env, stdlib))
    command = task.command.eval(env, stdlib).value.replace("/usr/local/src/", str(repo / "wdl/merge_results") + "/")
    subprocess.run(["bash", "-c", command], cwd=root, check=True, capture_output=True)
    assert (root / "rsem_isoforms_count.txt").read_text() == "transcript_id\tb\ta\ntx.1\t2.30\t1.20\n"
    assert all(path.is_symlink() for path in (root / "rsem_files").iterdir())
    assert task.runtime["memory"].eval(env, stdlib).value == "4GB"
    assert env.resolve("effective_scratch_gb").value == 11
print("isoform rendered command, sample order, symlinks, resource floors PASS")


task = WDL.load(str(repo / "wdl/rsem_exp/rsem.wdl")).tasks[0]
with tempfile.TemporaryDirectory() as directory:
    bam = Path(directory) / "input.bam"
    for gib, memory_floor, disk_floor, expected in (
        (0, 40, 60, (40, 60)), (0.06, 32, 60, (32, 60)), (12, 40, 60, (40, 60)),
        (15.21, 40, 60, (48, 71)), (30, 40, 60, (76, 130)),
        (15.21, 96, 200, (96, 200)), (70, 128, 250, (156, 290)),
    ):
        with bam.open("wb") as handle:
            handle.truncate(round(gib * 2**30))  # Sparse file exercises WDL size() without allocating data.
        env = WDL.values_from_json({"transcriptome_bam": str(bam), "memory": memory_floor,
                                    "disk_space": disk_floor}, task.available_inputs)
        for decl in task.postinputs:
            env = env.bind(decl.name, decl.expr.eval(env, Files("1.0")))
        assert task.runtime["memory"].eval(env, Files("1.0")).value == f"{expected[0]}GB"
        assert task.runtime["disks"].eval(env, Files("1.0")).value == f"local-disk {expected[1]} HDD"
print("RSEM growth and configured floors PASS")

workflow = WDL.load(str(repo / "wdl/rnaseq_pipeline_scatter.wdl")).workflow
scatter = next(node for node in workflow.body if isinstance(node, WDL.Tree.Scatter))
resources = [node for node in scatter.body if isinstance(node, WDL.Tree.Decl) and node.name in (
    "inferred_star_scratch_gb", "reference_star_scratch_gb", "effective_star_scratch_gb")]
for version, expected_tiers in (
    ("gencode_v47", [90, 120, 150, 180, 200, 250, 300, 400]),
    ("gencode_v50", [117, 156, 195, 234, 260, 325, 390, 520]),
    ("rn8", [90, 120, 150, 180, 200, 250, 300, 400]),
):
    for low, high, expected in zip(
        [0, 5000001, 40000001, 65000001, 90000001, 110000001, 155000001, 200000001],
        [5000000, 40000000, 65000000, 90000000, 110000000, 155000000, 200000000, 300000000],
        expected_tiers,
    ):
        for pairs in (low, high):
            for floor in (120, 600):
                env = WDL.Env.Bindings()
                for key, value in (("reference_release", WDL.Value.String(version)),
                                   ("cutadapt_read_pairs", WDL.Value.Int(pairs)),
                                   ("star_disk", WDL.Value.Int(floor))):
                    env = env.bind(key, value)
                for decl in resources:
                    env = env.bind(decl.name, decl.expr.eval(env, Files("1.0")))
                assert env.resolve("effective_star_scratch_gb").value == max(floor, expected)
defaults = {decl.name: decl.expr.eval(WDL.Env.Bindings(), Files("1.0")).value
            for decl in workflow.inputs if decl.name in (
                "num_preemptible_attempts", "use_umi_molecule_expression", "retain_all_read_expression")}
assert defaults == {"num_preemptible_attempts": 1, "use_umi_molecule_expression": True,
                    "retain_all_read_expression": False}
print("STAR v47/v50/rat tier boundaries, explicit floors, Spot and expression defaults PASS")

task = WDL.load(str(repo / "wdl/collect_rnaseq_metrics/collect_rnaseq_metrics.wdl")).tasks[0]
for enabled in (True, False):
    for cpu in (1, 2, 4):
        for ram in (8, 12, 13, 14, 48):
            env = WDL.values_from_json({"prefer_predefined_n1": enabled, "ncpu": cpu, "memory": ram}, task.available_inputs)
            selected = task.runtime["gcp"].eval(env, Files("1.0")).json
            expected = {"predefinedMachineType": "n1-highmem-2"} if enabled and cpu == 2 and ram == 12 else {}
            assert selected == expected
            assert task.runtime["cpu"].eval(env, Files("1.0")).value == cpu
            assert task.runtime["memory"].eval(env, Files("1.0")).value == f"{ram}GB"
print("Predefined N1 exact-shape selection, opt-out and larger resource requests PASS")
