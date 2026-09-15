"""Run the rendered UMI task on known paired-read families (MiniWDL + pysam)."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import tempfile

import pysam
import WDL


REPO = Path(__file__).resolve().parents[1]
# Three identical UMIs absorb one one-base error; a distant UMI stays separate.
FAMILIES = {"a0:AAAAAAAA", "a1:AAAAAAAA", "a2:AAAAAAAA"}
TEMPLATES = [(name, 100, 1) for name in sorted(FAMILIES)] + [
    ("error:AAAAAAAT", 100, 1), ("distinct:CCCCCCCC", 100, 1),
    ("multi:GGGGGGGG", 500, 2), ("invalid:NNNNNNNN", 700, 1),
]


class Files(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        return filename


def pair(name, start, nh, reference=0, secondary=False):
    for mate in (0, 1):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * 50
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.flag = (99 if mate == 0 else 147) | (256 if secondary else 0)
        read.reference_id = read.next_reference_id = reference
        read.reference_start = start + 100 * mate
        read.next_reference_start = start + 100 * (1 - mate)
        read.template_length = 150 if mate == 0 else -150
        read.mapping_quality = 60 if nh == 1 else 3
        read.cigarstring = "50M"
        read.set_tag("NH", nh)
        yield read


def make_bams(root):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 5000}]}
    records = [read for name, start, nh in TEMPLATES for read in pair(name, start, nh)]
    records += list(pair("multi:GGGGGGGG", 1000, 2, secondary=True))
    with pysam.AlignmentFile(str(root / "genomic.bam"), "wb", header=header) as bam:
        for read in sorted(records, key=lambda r: r.reference_start):
            bam.write(read)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": name, "LN": 5000} for name in ("tx1", "tx2")]}
    with pysam.AlignmentFile(str(root / "transcriptome.bam"), "wb", header=header) as bam:
        for name, _, _ in TEMPLATES:
            for reference in (0, 1):
                for read in pair(name, 50, 2, reference, secondary=reference == 1):
                    bam.write(read)


def read_records(path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return list(bam)


def check_outputs(root):
    genomic = read_records(root / "truth.umi_molecules.genomic.bam")
    selected = {read.query_name for read in genomic}
    assert len(selected & FAMILIES) == 1, selected
    assert selected - FAMILIES == {"distinct:CCCCCCCC", "multi:GGGGGGGG"}, selected
    for source, output in (("genomic.bam", "truth.umi_molecules.genomic.bam"),
                           ("transcriptome.bam", "truth.umi_molecules.transcriptome.bam")):
        expected = [r.to_string() for r in read_records(root / source) if r.query_name in selected]
        actual = [r.to_string() for r in read_records(root / output)]
        assert Counter(actual) == Counter(expected), output
    metrics = json.loads((root / "truth.umi_molecule_expression_metrics.json").read_text())
    denominators = metrics["denominators"]
    assert denominators["umi_eligible_genomic_templates"] == 6
    assert denominators["selected_genomic_representative_qnames"] == 3
    assert denominators["selected_genomic_alignment_records_for_featurecounts"] == 8
    assert denominators["selected_transcriptome_alignment_records_for_rsem"] == 12
    assert not list(root.glob("*.sqlite3*"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apptainer-image", type=Path, help="Local SIF of the pinned UMI base image")
    args = parser.parse_args()
    image = (REPO / "dockerfiles/umi_dup.Dockerfile").read_text().splitlines()[0].split()[1]
    task = WDL.load(str(REPO / "wdl/umi_dup/umi_dup.wdl")).tasks[0]
    with tempfile.TemporaryDirectory(prefix="umi-truth-") as directory:
        root = Path(directory).resolve()
        make_bams(root)
        env = WDL.values_from_json({"sample_prefix": "truth", "star_align": str(root / "genomic.bam"),
            "transcriptome_align": [str(root / "transcriptome.bam")], "emit_molecule_expression": True,
            "ncpu": 1, "memory": 2, "disk_space": 1, "preemptible": 0, "docker": image}, task.available_inputs)
        stdlib = Files("1.0", write_dir=directory)
        for decl in task.postinputs:
            env = env.bind(decl.name, decl.expr.eval(env, stdlib))
        command = task.command.eval(env, stdlib).value.replace("/usr/local/src/", str(REPO / "wdl/umi_dup") + "/")
        if args.apptainer_image:
            runner = ["apptainer", "exec", "--cleanenv", "--bind", f"{REPO}:{REPO}:ro",
                      "--bind", f"{root}:{root}", str(args.apptainer_image.resolve())]
        else:
            runner = ["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{REPO}:{REPO}:ro",
                      "-v", f"{root}:{root}", "-w", str(root), image]
        subprocess.run(runner + ["bash", "-c", command], cwd=root, check=True)
        check_outputs(root)
    print("UMI directional families, multimappers, transcript alternatives and denominators PASS")


if __name__ == "__main__":
    main()
