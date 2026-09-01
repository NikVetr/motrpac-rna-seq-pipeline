"""Merge featureCounts files in an explicit sample order."""

import argparse
import csv
from pathlib import Path


SUFFIX = ".out"


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
            raise ValueError("duplicate or invalid featureCounts sample: {}".format(sample))
        files[sample] = path
    if set(files) != set(samples):
        raise ValueError(
            "featureCounts samples differ from declared order: files={} order={}".format(
                sorted(files), sorted(samples)
            )
        )
    return files


def read_counts(path):
    with path.open(encoding="utf-8", newline="") as handle:
        first = handle.readline()
        if not first.startswith("#"):
            raise ValueError("featureCounts file lacks its metadata line: {}".format(path))
        rows = csv.reader(handle, delimiter="\t", strict=True)
        header = next(rows, None)
        if header is None or len(header) < 7 or header[0] != "Geneid":
            raise ValueError("invalid featureCounts header: {}".format(path))
        result = {}
        order = []
        for row in rows:
            if len(row) != len(header) or not row[0] or row[0] in result:
                raise ValueError("invalid featureCounts gene row: {}".format(path))
            order.append(row[0])
            result[row[0]] = row[6]
    if not order:
        raise ValueError("featureCounts file has no gene rows: {}".format(path))
    return order, result


def merge(directory, order_path, output_path):
    samples = sample_order(order_path)
    files = indexed_files(directory, samples)
    gene_order = None
    values = {}
    for sample in samples:
        current_order, current_values = read_counts(files[sample])
        if gene_order is None:
            gene_order = current_order
        if set(current_values) != set(gene_order):
            raise ValueError("featureCounts gene sets differ for {}".format(sample))
        values[sample] = current_values

    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id"] + samples)
        for gene in gene_order:
            writer.writerow([gene] + [values[sample][gene] for sample in samples])


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fc-dir", type=Path, required=True)
    parser.add_argument("--sample-order", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("featureCounts.txt"))
    args = parser.parse_args(argv)
    merge(args.fc_dir, args.sample_order, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
