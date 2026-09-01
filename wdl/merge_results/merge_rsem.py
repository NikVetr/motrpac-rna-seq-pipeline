"""Merge RSEM gene results in an explicit sample order."""

import argparse
import csv
from pathlib import Path


SUFFIX = ".genes.results"
METRICS = {
    "expected_count": "rsem_genes_count.txt",
    "TPM": "rsem_genes_tpm.txt",
    "FPKM": "rsem_genes_fpkm.txt",
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
            raise ValueError("duplicate or invalid RSEM sample: {}".format(sample))
        files[sample] = path
    if set(files) != set(samples):
        raise ValueError(
            "RSEM samples differ from declared order: files={} order={}".format(
                sorted(files), sorted(samples)
            )
        )
    return files


def read_results(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id"}.union(METRICS)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("RSEM file lacks required columns: {}".format(path))
        order = []
        values = {metric: {} for metric in METRICS}
        for row in reader:
            if None in row:
                raise ValueError("RSEM row has extra columns: {}".format(path))
            gene = row["gene_id"]
            if not gene or gene in values["expected_count"]:
                raise ValueError("invalid RSEM gene row: {}".format(path))
            order.append(gene)
            for metric in METRICS:
                if row[metric] in (None, ""):
                    raise ValueError("RSEM {} is missing for {}".format(metric, gene))
                values[metric][gene] = row[metric]
    if not order:
        raise ValueError("RSEM file has no gene rows: {}".format(path))
    return order, values


def merge(directory, order_path, output_directory):
    samples = sample_order(order_path)
    files = indexed_files(directory, samples)
    gene_order = None
    values = {}
    for sample in samples:
        current_order, current_values = read_results(files[sample])
        if gene_order is None:
            gene_order = current_order
        if set(current_values["expected_count"]) != set(gene_order):
            raise ValueError("RSEM gene sets differ for {}".format(sample))
        values[sample] = current_values

    for metric, filename in METRICS.items():
        output = output_directory / filename
        with output.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"] + samples)
            for gene in gene_order:
                writer.writerow(
                    [gene] + [values[sample][metric][gene] for sample in samples]
                )


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--rsem-dir", type=Path, required=True)
    parser.add_argument("--sample-order", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    merge(args.rsem_dir, args.sample_order, args.output_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
