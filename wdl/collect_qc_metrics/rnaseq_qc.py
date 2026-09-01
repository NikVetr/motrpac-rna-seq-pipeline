"""Collect the published per-sample QC row directly from native task reports."""

import argparse
import csv
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
import io
from pathlib import Path
import re
import tarfile
import zipfile


BASE_COLUMNS = [
    "sample",
    "reads_raw",
    "pct_adapter_detected",
    "pct_trimmed",
    "pct_trimmed_bases",
    "reads",
    "pct_GC",
    "pct_dup_sequence",
    "pct_rRNA",
    "pct_globin",
    "pct_phix",
    "pct_picard_dup",
]
STAR_COLUMNS = [
    "avg_input_read_length",
    "uniquely_mapped",
    "pct_uniquely_mapped",
    "avg_mapped_read_length",
    "num_splices",
    "num_annotated_splices",
    "num_GTAG_splices",
    "num_GCAG_splices",
    "num_ATAC_splices",
    "num_noncanonical_splices",
    "pct_multimapped",
    "pct_multimapped_toomany",
    "pct_unmapped_mismatches",
    "pct_unmapped_tooshort",
    "pct_unmapped_other",
    "pct_chimeric",
]
MAPPED_COLUMNS = ["pct_chrX", "pct_chrY", "pct_chrM", "pct_chrAuto", "pct_contig"]
RNA_COLUMNS = [
    "pct_coding",
    "pct_utr",
    "pct_intronic",
    "pct_intergenic",
    "pct_mrna",
    "median_5_3_bias",
]


def read_text(path):
    text = path.read_text(encoding="utf-8")
    if not text:
        raise ValueError("input report is empty: {}".format(path))
    return text


def number(value, label, percent=False):
    text = str(value).strip().replace(",", "")
    if text.endswith("%"):
        text = text[:-1]
    try:
        result = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError("{} is not numeric: {!r}".format(label, value)) from exc
    if not result.is_finite() or result < 0:
        raise ValueError("{} must be finite and nonnegative".format(label))
    if percent and result > 100:
        raise ValueError("{} is outside 0-100".format(label))
    return result


def integer(value, label):
    result = number(value, label)
    if result != result.to_integral_value():
        raise ValueError("{} is not an integer".format(label))
    return int(result)


def rendered(value, places=None):
    if not isinstance(value, Decimal):
        value = Decimal(value)
    if places is not None:
        value = value.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_EVEN)
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text if text and text != "-0" else "0"


def single_match(pattern, text, label):
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise ValueError("expected one {} value; found {}".format(label, len(matches)))
    return matches[0]


def parse_fastqc_data(data, source):
    text = data.decode("utf-8")
    fields = {}
    for line in text.splitlines():
        if line.startswith(("Filename\t", "Total Sequences\t", "%GC\t")):
            key, value = line.split("\t", 1)
            if key in fields:
                raise ValueError("duplicate FastQC {} in {}".format(key, source))
            fields[key] = value
    for key in ("Filename", "Total Sequences", "%GC"):
        if key not in fields:
            raise ValueError("FastQC {} lacks {}".format(source, key))
    deduplicated = number(
        single_match(
            r"^#Total Deduplicated Percentage\t([^\t]+)$",
            text,
            "FastQC deduplicated percentage",
        ),
        "FastQC deduplicated percentage",
        percent=True,
    )
    return {
        "filename": fields["Filename"],
        "reads": integer(fields["Total Sequences"], "FastQC total sequences"),
        "gc": number(fields["%GC"], "FastQC GC", percent=True),
        "duplicates": Decimal(100) - deduplicated,
    }


def parse_fastqc_archive(path, expected_filenames):
    reports = {}
    with tarfile.open(str(path), "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.endswith("_fastqc.zip")
        ]
        if len(members) != 2:
            raise ValueError("{} must contain two FastQC ZIP reports".format(path))
        for member in members:
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError("could not read {} from {}".format(member.name, path))
            with zipfile.ZipFile(io.BytesIO(extracted.read())) as zipped:
                candidates = [
                    name for name in zipped.namelist() if name.endswith("/fastqc_data.txt")
                ]
                if len(candidates) != 1:
                    raise ValueError("{} lacks one fastqc_data.txt".format(member.name))
                report = parse_fastqc_data(zipped.read(candidates[0]), member.name)
            filename = report["filename"]
            if filename in reports:
                raise ValueError("duplicate FastQC filename: {}".format(filename))
            reports[filename] = report
    if set(reports) != set(expected_filenames):
        raise ValueError(
            "FastQC filenames {} differ from expected {}".format(
                sorted(reports), sorted(expected_filenames)
            )
        )
    return reports


def parse_cutadapt(path):
    text = read_text(path)
    patterns = {
        "pairs_in": r"^Total read pairs processed:\s+([0-9,]+)\s*$",
        "r1_adapter_pct": r"^\s*Read 1 with adapter:\s+[0-9,]+\s+\(([0-9.]+)%\)\s*$",
        "r2_adapter_pct": r"^\s*Read 2 with adapter:\s+[0-9,]+\s+\(([0-9.]+)%\)\s*$",
        "pairs_out": r"^Pairs written \(passing filters\):\s+([0-9,]+)\s+\([^\n]+\)\s*$",
        "bases_in": r"^Total basepairs processed:\s+([0-9,]+)\s+bp\s*$",
        "bases_out": r"^Total written \(filtered\):\s+([0-9,]+)\s+bp(?:\s+\([^\n]+\))?\s*$",
    }
    result = {}
    for key, pattern in patterns.items():
        raw = single_match(pattern, text, "Cutadapt {}".format(key))
        result[key] = (
            number(raw, "Cutadapt {}".format(key), percent=True)
            if key.endswith("_pct")
            else integer(raw, "Cutadapt {}".format(key))
        )
    if result["pairs_out"] > result["pairs_in"] or result["bases_out"] > result["bases_in"]:
        raise ValueError("Cutadapt output counts exceed input counts")
    return result


def single_row(path, delimiter="\t"):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
    if reader.fieldnames is None or len(rows) != 1 or None in rows[0]:
        raise ValueError("{} must contain one complete data row".format(path))
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("{} has duplicate columns".format(path))
    if any(rows[0].get(name) in (None, "") for name in reader.fieldnames):
        raise ValueError("{} has missing values".format(path))
    return rows[0]


def sample_percentage(path, sample, metric):
    row = single_row(path)
    if row.get("Sample") != sample or metric not in row:
        raise ValueError("{} does not identify {} with {}".format(path, sample, metric))
    number(row[metric], metric, percent=True)
    return row[metric].strip().rstrip("%")


def parse_star(path):
    values = {}
    for line in read_text(path).splitlines():
        if "|" not in line:
            continue
        key, value = (part.strip() for part in line.split("|", 1))
        if not key or key in values:
            raise ValueError("STAR log has an empty or duplicate field")
        values[key] = value
    return values


def required(mapping, key, source):
    if mapping.get(key, "") == "":
        raise ValueError("{} lacks {}".format(source, key))
    return mapping[key]


def parse_picard(path, metrics_class):
    lines = read_text(path).splitlines()
    marker = "## METRICS CLASS\t{}".format(metrics_class)
    positions = [index for index, line in enumerate(lines) if line == marker]
    if len(positions) != 1:
        raise ValueError("{} lacks one {} section".format(path, metrics_class))
    records = []
    for line in lines[positions[0] + 1 :]:
        if not line:
            if records:
                break
            continue
        if line.startswith("#"):
            if records:
                break
            continue
        records.append(line)
    if len(records) != 2:
        raise ValueError("{} must contain one {} row".format(path, metrics_class))
    header, values = records[0].split("\t"), records[1].split("\t")
    if len(header) != len(set(header)) or len(values) != len(header):
        raise ValueError("malformed Picard metrics: {}".format(path))
    return dict(zip(header, values))


def fraction_percent(metrics, key, source):
    value = number(required(metrics, key, source), "{} {}".format(source, key))
    if value > 1:
        raise ValueError("{} {} is outside 0-1".format(source, key))
    return value * 100


def collect(args):
    cutadapt = parse_cutadapt(args.cutadapt_report)
    reads_raw, reads = cutadapt["pairs_in"], cutadapt["pairs_out"]
    if reads_raw == 0:
        raise ValueError("raw read-pair count is zero")
    if cutadapt["bases_in"] == 0:
        raise ValueError("Cutadapt input base count is zero")

    pre_reports = None
    if args.fastqc_pretrim is not None:
        pre = parse_fastqc_archive(
            args.fastqc_pretrim, (args.pretrim_r1_filename, args.pretrim_r2_filename)
        )
        pre_reports = [pre[args.pretrim_r1_filename], pre[args.pretrim_r2_filename]]
        if pre_reports[0]["reads"] != pre_reports[1]["reads"]:
            raise ValueError("pre-trim FastQC mate counts disagree")
        if pre_reports[0]["reads"] != reads_raw:
            raise ValueError("Cutadapt and pre-trim FastQC read counts disagree")

    post_reports = None
    if args.fastqc_posttrim is not None:
        post = parse_fastqc_archive(
            args.fastqc_posttrim,
            (args.posttrim_r1_filename, args.posttrim_r2_filename),
        )
        post_reports = [post[args.posttrim_r1_filename], post[args.posttrim_r2_filename]]
        if post_reports[0]["reads"] != post_reports[1]["reads"]:
            raise ValueError("post-trim FastQC mate counts disagree")
        if post_reports[0]["reads"] != reads:
            raise ValueError("Cutadapt and post-trim FastQC read counts disagree")

    star = parse_star(args.star_log)
    if integer(required(star, "Number of input reads", "STAR"), "STAR input reads") != reads:
        raise ValueError("STAR and Cutadapt read counts disagree")
    markdup = (
        parse_picard(args.markduplicates_metrics, "picard.sam.DuplicationMetrics")
        if args.markduplicates_metrics is not None
        else None
    )
    rna = (
        parse_picard(args.rnaseq_metrics, "picard.analysis.RnaSeqMetrics")
        if args.rnaseq_metrics is not None
        else None
    )
    mapped = single_row(args.mapped_report) if args.mapped_report is not None else None
    if mapped is not None and mapped.get("Sample") != args.sample:
        raise ValueError("mapped report sample differs from authoritative sample")

    row = {
        "sample": args.sample,
        "reads_raw": str(reads_raw),
        "pct_adapter_detected": rendered(
            (cutadapt["r1_adapter_pct"] + cutadapt["r2_adapter_pct"]) / 2, 3
        ),
        "pct_trimmed": rendered(Decimal(reads_raw - reads) / Decimal(reads_raw) * 100, 3),
        "pct_trimmed_bases": rendered(
            Decimal(cutadapt["bases_in"] - cutadapt["bases_out"])
            / Decimal(cutadapt["bases_in"])
            * 100,
            3,
        ),
        "reads": str(reads),
        "pct_GC": (
            rendered((post_reports[0]["gc"] + post_reports[1]["gc"]) / 2, 3)
            if post_reports is not None
            else ""
        ),
        "pct_dup_sequence": (
            rendered(
                (post_reports[0]["duplicates"] + post_reports[1]["duplicates"]) / 2,
                3,
            )
            if post_reports is not None
            else ""
        ),
        "pct_rRNA": (
            sample_percentage(args.rrna_report, args.sample, "pct_rRNA")
            if args.rrna_report is not None
            else ""
        ),
        "pct_globin": (
            sample_percentage(args.globin_report, args.sample, "pct_globin")
            if args.globin_report is not None
            else ""
        ),
        "pct_phix": (
            sample_percentage(args.phix_report, args.sample, "pct_phix")
            if args.phix_report is not None
            else ""
        ),
        "pct_picard_dup": (
            rendered(
                fraction_percent(
                    markdup,
                    "PERCENT_DUPLICATION",
                    str(args.markduplicates_metrics),
                ),
                3,
            )
            if markdup is not None
            else ""
        ),
    }
    row["pct_umi_dup"] = (
        sample_percentage(args.umi_report, args.sample, "pct_umi_dup")
        if args.umi_report is not None
        else ""
    )

    star_fields = {
        "avg_input_read_length": ("Average input read length", False, False),
        "uniquely_mapped": ("Uniquely mapped reads number", False, True),
        "pct_uniquely_mapped": ("Uniquely mapped reads %", True, False),
        "avg_mapped_read_length": ("Average mapped length", False, False),
        "num_splices": ("Number of splices: Total", False, True),
        "num_annotated_splices": ("Number of splices: Annotated (sjdb)", False, True),
        "num_GTAG_splices": ("Number of splices: GT/AG", False, True),
        "num_GCAG_splices": ("Number of splices: GC/AG", False, True),
        "num_ATAC_splices": ("Number of splices: AT/AC", False, True),
        "num_noncanonical_splices": ("Number of splices: Non-canonical", False, True),
        "pct_multimapped": ("% of reads mapped to multiple loci", True, False),
        "pct_multimapped_toomany": ("% of reads mapped to too many loci", True, False),
        "pct_unmapped_mismatches": ("% of reads unmapped: too many mismatches", True, False),
        "pct_unmapped_tooshort": ("% of reads unmapped: too short", True, False),
        "pct_unmapped_other": ("% of reads unmapped: other", True, False),
    }
    for output_name, (source_name, percent, integral) in star_fields.items():
        raw = required(star, source_name, "STAR")
        if integral:
            integer(raw, "STAR {}".format(source_name))
        else:
            number(raw, "STAR {}".format(source_name), percent=percent)
        row[output_name] = raw.rstrip("%")
    # STAR chimeric detection is disabled, so its reported zero is unavailable data.
    row["pct_chimeric"] = ""

    mapped_places = {"pct_chrX": 3, "pct_chrY": 5, "pct_chrM": 3, "pct_chrAuto": 3, "pct_contig": 3}
    for name, places in mapped_places.items():
        row[name] = (
            rendered(number(required(mapped, name, "mapped"), name, percent=True), places)
            if mapped is not None
            else ""
        )

    rna_fields = {
        "pct_coding": "PCT_CODING_BASES",
        "pct_utr": "PCT_UTR_BASES",
        "pct_intronic": "PCT_INTRONIC_BASES",
        "pct_intergenic": "PCT_INTERGENIC_BASES",
        "pct_mrna": "PCT_MRNA_BASES",
    }
    for output_name, source_name in rna_fields.items():
        row[output_name] = (
            rendered(fraction_percent(rna, source_name, str(args.rnaseq_metrics)), 3)
            if rna is not None
            else ""
        )
    row["median_5_3_bias"] = (
        rendered(
            number(
                required(rna, "MEDIAN_5PRIME_TO_3PRIME_BIAS", "RNA metrics"),
                "median bias",
            ),
            3,
        )
        if rna is not None
        else ""
    )

    columns = BASE_COLUMNS + ["pct_umi_dup"] + STAR_COLUMNS + MAPPED_COLUMNS + RNA_COLUMNS
    if list(row) != columns:
        raise AssertionError("QC fields do not match the published column contract")
    return columns, row


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fastqc-pretrim", type=Path)
    parser.add_argument("--fastqc-posttrim", type=Path)
    parser.add_argument("--pretrim-r1-filename", required=True)
    parser.add_argument("--pretrim-r2-filename", required=True)
    parser.add_argument("--posttrim-r1-filename", required=True)
    parser.add_argument("--posttrim-r2-filename", required=True)
    parser.add_argument("--cutadapt-report", type=Path, required=True)
    parser.add_argument("--mapped-report", type=Path)
    parser.add_argument("--rrna-report", type=Path)
    parser.add_argument("--globin-report", type=Path)
    parser.add_argument("--phix-report", type=Path)
    parser.add_argument("--star-log", type=Path, required=True)
    parser.add_argument("--markduplicates-metrics", type=Path)
    parser.add_argument("--rnaseq-metrics", type=Path)
    parser.add_argument("--umi-report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv=None):
    args = make_parser().parse_args(argv)
    columns, row = collect(args)
    with args.output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
