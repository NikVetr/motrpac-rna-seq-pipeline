#!/usr/bin/env python3
"""Select the fixed production UMI grouping population and attach RX tags."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path


UMI_LENGTH = 8
UMI_PATTERN = re.compile(r"^[ACGTN]{8}$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bam", type=Path, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    return parser.parse_args()


def umi_from_query_name(query_name: str) -> str:
    if ":" not in query_name:
        raise ValueError(f"read name lacks colon-delimited UMI: {query_name}")
    umi = query_name.rsplit(":", 1)[1].upper()
    if not UMI_PATTERN.fullmatch(umi):
        raise ValueError(
            f"read name has invalid {UMI_LENGTH}-base UMI {umi!r}: {query_name}"
        )
    return umi


def main() -> int:
    import pysam

    args = parse_args()
    if not args.input_bam.is_file():
        raise FileNotFoundError(f"input BAM does not exist: {args.input_bam}")

    counters: Counter[str] = Counter()
    with pysam.AlignmentFile(str(args.input_bam), "rb") as source:
        with pysam.AlignmentFile(str(args.output_bam), "wb", template=source) as output:
            for read in source.fetch(until_eof=True):
                counters["input_alignment_records"] += 1
                umi = umi_from_query_name(read.query_name)

                if read.is_secondary:
                    counters["excluded_secondary_records"] += 1
                    continue
                if read.is_supplementary:
                    counters["excluded_supplementary_records"] += 1
                    continue
                if not read.is_paired:
                    counters["excluded_unpaired_records"] += 1
                    continue
                if read.is_unmapped or read.mate_is_unmapped:
                    counters["excluded_unmapped_records"] += 1
                    continue
                if not read.is_proper_pair:
                    counters["excluded_chimeric_or_improper_records"] += 1
                    continue
                if not (read.is_read1 or read.is_read2):
                    raise ValueError(
                        f"eligible paired record lacks an R1/R2 flag: {read.query_name}"
                    )
                if not read.has_tag("NH"):
                    raise ValueError(
                        f"eligible record lacks STAR NH tag: {read.query_name}"
                    )

                mate = "r1" if read.is_read1 else "r2"
                counters[f"primary_proper_{mate}_records"] += 1
                if "N" in umi:
                    counters[f"excluded_n_umi_{mate}_records"] += 1
                    continue

                nh = read.get_tag("NH")
                if isinstance(nh, bool) or not isinstance(nh, int) or nh < 1:
                    raise ValueError(
                        f"eligible record has invalid STAR NH tag {nh!r}: "
                        f"{read.query_name}"
                    )
                mapping_class = "unique" if nh == 1 else "multimapped"
                counters[f"eligible_{mapping_class}_{mate}_records"] += 1
                counters[f"eligible_acgt_{mate}_records"] += 1
                counters["eligible_alignment_records"] += 1
                read.set_tag("RX", umi, value_type="Z", replace=True)
                output.write(read)

    for category in (
        "primary_proper",
        "excluded_n_umi",
        "eligible_unique",
        "eligible_multimapped",
        "eligible_acgt",
    ):
        r1 = counters[f"{category}_r1_records"]
        r2 = counters[f"{category}_r2_records"]
        if r1 != r2:
            raise ValueError(f"{category} R1/R2 record counts differ: {r1} != {r2}")

    eligible = counters["eligible_acgt_r1_records"]
    if eligible == 0:
        raise ValueError("no eligible A/C/G/T UMI templates remain")
    if counters["eligible_alignment_records"] != 2 * eligible:
        raise ValueError("eligible record and template counts disagree")

    primary = counters["primary_proper_r1_records"]
    invalid_n = counters["excluded_n_umi_r1_records"]
    nh1 = counters["eligible_unique_r1_records"]
    nh_gt1 = counters["eligible_multimapped_r1_records"]
    if primary != invalid_n + nh1 + nh_gt1:
        raise ValueError("primary-proper UMI denominator partition is inconsistent")

    metrics = {
        "schema_version": 1,
        "representation": "rx_v1",
        "umi_length": UMI_LENGTH,
        "mapping_eligibility": "all_primary_proper",
        "policy": {
            "alignment_records": (
                "primary proper pairs only; secondary, supplementary, unmapped, "
                "unpaired, and improper/chimeric records excluded"
            ),
            "mapping_quality": (
                "STAR NH tag required; all primary proper A/C/G/T templates eligible"
            ),
            "umi": "exactly eight bases; A/C/G/T templates eligible; N templates excluded",
        },
        **dict(sorted(counters.items())),
        "primary_proper_templates": primary,
        "excluded_n_umi_templates": invalid_n,
        "eligible_acgt_templates": eligible,
        "eligible_nh1_templates": nh1,
        "eligible_multimapped_templates": nh_gt1,
    }
    args.metrics.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
