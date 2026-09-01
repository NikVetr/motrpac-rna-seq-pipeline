"""Gather one-row QC reports in an explicit sample order."""

import argparse
import csv
from pathlib import Path


SUFFIX = "_qc_info.csv"
OPTIONAL_QC_COLUMNS = {
    "pct_GC",
    "pct_dup_sequence",
    "pct_rRNA",
    "pct_globin",
    "pct_phix",
    "pct_picard_dup",
    "pct_umi_dup",
    "pct_chimeric",
    "pct_chrX",
    "pct_chrY",
    "pct_chrM",
    "pct_chrAuto",
    "pct_contig",
    "pct_coding",
    "pct_utr",
    "pct_intronic",
    "pct_intergenic",
    "pct_mrna",
    "median_5_3_bias",
}


def sample_order(path):
    samples = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if not samples or any(not sample for sample in samples):
        raise ValueError("sample order must contain nonempty sample IDs")
    if len(samples) != len(set(samples)):
        raise ValueError("sample order contains duplicate IDs")
    return samples


def indexed_files(directory, samples):
    files = {}
    for path in directory.iterdir():
        if not path.is_file() or not path.name.endswith(SUFFIX):
            continue
        sample = path.name[: -len(SUFFIX)]
        if not sample or sample in files:
            raise ValueError("duplicate or invalid QC sample: {}".format(sample))
        files[sample] = path
    if set(files) != set(samples):
        raise ValueError(
            "QC samples differ from declared order: files={} order={}".format(
                sorted(files), sorted(samples)
            )
        )
    return files


def read_row(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if reader.fieldnames is None or len(rows) != 1:
        raise ValueError("QC report must contain one header and one row: {}".format(path))
    if len(reader.fieldnames) != len(set(reader.fieldnames)) or "sample" not in reader.fieldnames:
        raise ValueError("QC report has invalid columns: {}".format(path))
    if None in rows[0] or any(
        rows[0].get(name) in (None, "")
        for name in reader.fieldnames
        if name not in OPTIONAL_QC_COLUMNS
    ):
        raise ValueError("QC report has missing values: {}".format(path))
    return reader.fieldnames, rows[0]


def merge(directory, order_path, output_path):
    samples = sample_order(order_path)
    files = indexed_files(directory, samples)
    header = None
    rows = []
    for sample in samples:
        current_header, row = read_row(files[sample])
        if header is None:
            header = current_header
        if current_header != header:
            raise ValueError("QC columns differ for {}".format(sample))
        if row["sample"] != sample:
            raise ValueError("QC row sample differs from filename: {}".format(sample))
        rows.append(row)

    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--sample-order", type=Path, required=True)
    parser.add_argument("--output-name", type=Path, required=True)
    args = parser.parse_args(argv)
    merge(args.qc_dir, args.sample_order, args.output_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
