"""Join study covariates and modern QC in expression-column order without rescaling."""
import argparse
import csv
from pathlib import Path


def read_table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if not columns or "sample" not in columns or len(columns) != len(set(columns)):
            raise ValueError(f"invalid metadata header: {path}")
        rows = {}
        for row in reader:
            sample = row["sample"]
            if not sample or sample in rows or None in row or None in row.values():
                raise ValueError(f"duplicate sample or malformed metadata row: {path}")
            rows[sample] = row
    return columns, rows


def prepare(matrix, study, qc, output, required):
    with matrix.open(encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    samples = header[1:]
    if header[0] != "gene_id" or not samples or "" in samples or len(samples) != len(set(samples)):
        raise ValueError("invalid expression sample header")
    study_columns, study_rows = read_table(study)
    qc_columns, qc_rows = read_table(qc)
    if set(study_columns) & set(qc_columns) != {"sample"}:
        raise ValueError("study and QC columns overlap; remove historical pipeline QC from the study table")
    if not set(required).issubset(study_columns):
        raise ValueError("required study columns are missing")
    if not set(samples).issubset(study_rows) or set(samples) != set(qc_rows):
        raise ValueError("study or QC samples do not match expression columns")
    for sample in samples:
        if any(study_rows[sample][column].strip().upper() in ("", "NA", "NAN", "NULL") for column in required):
            raise ValueError(f"missing required study covariate: {sample}")
    with output.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=study_columns + [column for column in qc_columns if column != "sample"])
        writer.writeheader()
        writer.writerows({**study_rows[sample], **qc_rows[sample]} for sample in samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("matrix", "study", "qc", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--required-study-columns", nargs="+", required=True)
    args = parser.parse_args()
    prepare(args.matrix, args.study, args.qc, args.output, args.required_study_columns)


if __name__ == "__main__":
    main()
