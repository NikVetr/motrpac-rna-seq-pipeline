"""Merge RSEM gene or isoform results in an explicit sample order."""

import argparse
import csv
from contextlib import ExitStack
from itertools import zip_longest
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


def indexed_files(directory, samples, suffix=SUFFIX):
    files = {}
    for path in directory.iterdir():
        if not path.is_file() or not path.name.endswith(suffix):
            continue
        sample = path.name[: -len(suffix)]
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


def merge_isoforms(directory, order_path, output_directory):
    samples = sample_order(order_path)
    files = indexed_files(directory, samples, ".isoforms.results")
    with ExitStack() as stack:
        readers = [csv.DictReader(stack.enter_context(files[sample].open(
            encoding="utf-8", newline="")), delimiter="\t") for sample in samples]
        required = {"transcript_id", "gene_id"}.union(METRICS)
        for sample, reader in zip(samples, readers):
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError("RSEM isoform file lacks required columns: {}".format(sample))
        writers = {}
        for metric, filename in METRICS.items():
            handle = stack.enter_context((output_directory / filename.replace("genes", "isoforms")).open(
                "x", encoding="utf-8", newline=""))
            writers[metric] = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writers[metric].writerow(["transcript_id"] + samples)
        seen = set()
        for rows in zip_longest(*readers):
            if any(row is None or None in row or any(row.get(key) in (None, "") for key in required)
                   for row in rows):
                raise ValueError("RSEM isoform rows are malformed or have different lengths")
            transcript = rows[0]["transcript_id"]
            if transcript in seen:
                raise ValueError("duplicate RSEM transcript: {}".format(transcript))
            seen.add(transcript)
            if any((row["transcript_id"], row["gene_id"]) != (transcript, rows[0]["gene_id"]) for row in rows):
                raise ValueError("RSEM transcript order or gene mapping differs: {}".format(transcript))
            for metric, writer in writers.items():
                writer.writerow([transcript] + [row[metric] for row in rows])
        if not seen:
            raise ValueError("RSEM isoform files have no transcript rows")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--rsem-dir", type=Path, required=True)
    parser.add_argument("--sample-order", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, default=Path("."))
    parser.add_argument("--feature-level", choices=("genes", "isoforms"), default="genes")
    args = parser.parse_args(argv)
    merger = merge_isoforms if args.feature_level == "isoforms" else merge
    merger(args.rsem_dir, args.sample_order, args.output_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
