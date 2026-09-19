#!/usr/bin/env python3
"""Exercise the diagnostic writer against RSEM 1.3.3's actual transcript/group types."""
import argparse
import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = r'''
#include "rsem_convergence.h"
int main() {
    Transcripts transcripts;
    for (int i = 1; i <= 4; ++i)
        transcripts.add(Transcript("t" + std::to_string(i), i <= 2 ? "g1" : "g2",
            "chr1", '+', std::vector<Interval>{Interval(1, 100)}, ""));
    std::ofstream("reference.grp") << "1\n3\n5\n";
    std::vector<double> before{0.1, 0.0001, 0.4, 1e-8, 1e-7};
    std::vector<double> after{0.1002, 0.000101, 0.399999, 1e-7, 2e-7};
    std::vector<double> oldCounts{10, 0.01, 40, 0.000001, 0.1};
    std::vector<double> counts{10.5, 0.11, 39.9, 0.000002, 0.2};
    writeConvergence("capped", "reference", transcripts, 20, 0.001,
                     before, after, oldCounts, counts.data());
    writeConvergence("converged", "reference", transcripts, 21, 0.001,
                     after, after, counts, counts.data());
    try {
        writeConvergence("missing/directory/sample", "reference", transcripts, 20,
                         0.001, before, after, oldCounts, counts.data());
        return 1;
    } catch (const std::ios_base::failure&) {}
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rsem-source", type=Path, required=True)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("qc", ROOT / "wdl/collect_qc_metrics/rnaseq_qc.py")
    qc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(qc)
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        source = work / "fixture.cpp"
        source.write_text(FIXTURE)
        subprocess.run([args.compiler, "-std=c++14", "-O3", "-ffast-math",
                        "-I" + str(ROOT / "dockerfiles"), "-I" + str(args.rsem_source.resolve()),
                        str(source), "-o", str(work / "fixture")], check=True, capture_output=True)
        subprocess.run([str(work / "fixture")], cwd=work, check=True)
        for name, remaining in (("capped", 3), ("converged", 0)):
            result = qc.parse_rsem_convergence(work / (name + ".rsem_convergence.tsv"),
                work / (name + ".rsem_gene_convergence.tsv"),
                {"iterations": 20 if remaining else 21, "components_above_tolerance": remaining})
            assert result["flagged_transcripts"] == (2 if remaining else 0)
            assert result["affected_genes"] == (2 if remaining else 0)
            assert (result["background"] is not None) == bool(remaining)
        with (work / "capped.rsem_convergence.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["transcript_id"] for row in rows] == [".", "t1", "t4"]
        with (work / "capped.rsem_gene_convergence.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        assert genes[0]["flagged_transcripts"] == "1"
        assert math.isclose(float(genes[0]["expected_count_final_theta"]), 40.01)
        assert math.isclose(float(genes[0]["expected_count_change"]), 0, abs_tol=1e-12)
        assert math.isclose(float(genes[0]["sum_absolute_isoform_count_change"]), 0.2)
    print("PASS: thresholds, background, all-isoform gene sums, empty reports and write failures")


if __name__ == "__main__":
    main()
