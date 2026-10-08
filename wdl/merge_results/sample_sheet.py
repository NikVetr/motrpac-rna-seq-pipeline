"""Write the per-sample sheet: expression metadata plus status, diagnostics, input sizes and QC flags."""

import argparse
import csv
import json
from pathlib import Path


SIZE_COLUMNS = ["raw_fastq_gib", "genomic_bam_gib", "transcriptome_bam_gib", "rsem_input_bam_gib"]
COLUMNS = ["status", "rsem_converged", "rsem_iterations", "correct_strand_fraction",
           "featurecounts_assigned_fraction"] + SIZE_COLUMNS + ["qc_flags"]


def number(value):
    return None if value in (None, "") else float(value)


def qc_flags(qc, diagnostics, cohort_mean_mapped):
    """Review-only flags: GET MOP section 9 thresholds plus strandedness and RSEM convergence."""
    # Absent or blank operands (e.g. a skipped QC group) leave that flag unevaluated.
    reads, gc, rrna = number(qc.get("reads")), number(qc.get("pct_GC")), number(qc.get("pct_rRNA"))
    unique, multi = number(qc.get("pct_uniquely_mapped")), number(qc.get("pct_multimapped"))
    coding, utr = number(qc.get("pct_coding")), number(qc.get("pct_utr"))
    mapped = None if None in (unique, multi) else unique + multi
    strand = (diagnostics.get("picard_strand") or {}).get("correct_strand_fraction")
    checks = {
        "low_reads_after_trim": reads is not None and reads < 20e6,
        "abnormal_gc": gc is not None and not 20 <= gc <= 80,
        "high_rrna": rrna is not None and rrna > 20,
        "low_mapped": mapped is not None and mapped < 60,
        "low_exonic": None not in (coding, utr) and coding + utr < 50,
        "low_mapped_vs_cohort": None not in (reads, mapped, cohort_mean_mapped)
                                and reads * mapped / 100 < 0.5 * cohort_mean_mapped,
        "low_strand": strand is not None and strand < 0.9,
        "rsem_not_converged": not diagnostics["rsem"]["converged"],
    }
    return ";".join(name for name, flagged in checks.items() if flagged)


def read_tsv(path, header=True):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows if not header else (rows[0], rows[1:]) if rows else ([], [])


def build(metadata, qc_report, diagnostics_dir, sizes, failed, output):
    header, rows = read_tsv(metadata)
    samples = [row[0] for row in rows]
    failed = {line.strip() for line in failed.read_text(encoding="utf-8").splitlines() if line.strip()}
    if header[:1] != ["sample"] or len(samples) != len(set(samples)) or not failed <= set(samples):
        raise ValueError("metadata must list unique samples, including every failed sample")
    completed = [sample for sample in samples if sample not in failed]
    with qc_report.open(encoding="utf-8", newline="") as handle:
        qc = {row["sample"]: row for row in csv.DictReader(handle)}
    diagnostics = {}
    for path in diagnostics_dir.glob("*.qc_diagnostics.json"):
        record = json.loads(path.read_text(encoding="utf-8"))
        diagnostics[record["sample"]] = record
    size_rows = {row[0]: row[1:] for row in read_tsv(sizes, header=False) if row}
    if set(qc) != set(completed) or set(diagnostics) != set(completed):
        raise ValueError("QC report and diagnostics must cover exactly the completed samples")
    if size_rows and (set(size_rows) != set(completed) or
                      any(len(values) != len(SIZE_COLUMNS) for values in size_rows.values())):
        raise ValueError("size rows must cover exactly the completed samples")

    operands = [[number(qc[s].get(key)) for key in ("reads", "pct_uniquely_mapped", "pct_multimapped")] for s in completed]
    mapped = [reads * (unique + multi) / 100 for reads, unique, multi in operands if None not in (reads, unique, multi)]
    cohort_mean_mapped = sum(mapped) / len(mapped) if mapped else None
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header + COLUMNS)
        for row in rows:
            sample = row[0]
            if sample in failed:
                writer.writerow(row + ["failed"] + [""] * (len(COLUMNS) - 1))
                continue
            record = diagnostics[sample]
            strand = (record.get("picard_strand") or {}).get("correct_strand_fraction")
            sizes_text = ["" if value == "" else "{:.3f}".format(float(value))
                          for value in size_rows.get(sample, [""] * len(SIZE_COLUMNS))]
            writer.writerow(row + [
                "completed", str(record["rsem"]["converged"]).lower(), record["rsem"]["iterations"],
                "" if strand is None else "{:.6f}".format(strand),
                "{:.6f}".format(record["feature_counts"]["assigned_alignment_fraction"]),
                *sizes_text, qc_flags(qc[sample], record, cohort_mean_mapped)])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True, help="expression metadata rows with header")
    parser.add_argument("--qc", type=Path, required=True, help="consolidated QC CSV")
    parser.add_argument("--diagnostics-dir", type=Path, required=True)
    parser.add_argument("--sizes", type=Path, required=True, help="headerless rows: sample + SIZE_COLUMNS; may be empty")
    parser.add_argument("--failed", type=Path, required=True, help="failed sample IDs, one per line; may be empty")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    build(args.metadata, args.qc, args.diagnostics_dir, args.sizes, args.failed, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
