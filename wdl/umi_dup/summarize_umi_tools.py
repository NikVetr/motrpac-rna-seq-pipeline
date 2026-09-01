#!/usr/bin/env python3
"""Validate UMI-tools output and write production QC and structured metrics."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path


EXPECTED_VERSION = "UMI-tools version: 1.1.6"
PATTERNS = {
    "tool_input_templates": re.compile(r"Reads: Input Reads: (\d+), Read pairs: (\d+)"),
    "molecules": re.compile(r"Number of reads out: (\d+)"),
    "positions": re.compile(r"Total number of positions deduplicated: (\d+)"),
    "mean_unique_umis_per_position": re.compile(
        r"Mean number of unique UMIs per position: ([0-9.eE+-]+)"
    ),
    "max_unique_umis_per_position": re.compile(
        r"Max\. number of unique UMIs per position: (\d+)"
    ),
    "unmatched_alignments": re.compile(
        r"Searching for mates for (\d+) unmatched alignments"
    ),
    "mates_never_found": re.compile(r"(\d+) mates never found"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--container", required=True)
    return parser.parse_args()


def exactly_one(pattern: re.Pattern[str], text: str, name: str) -> tuple[str, ...]:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise ValueError(
            f"expected one {name} record in UMI-tools log, found {len(matches)}"
        )
    value = matches[0]
    return value if isinstance(value, tuple) else (value,)


def rendered_percentage(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def main() -> int:
    args = parse_args()
    if args.version.strip() != EXPECTED_VERSION:
        raise ValueError(f"unexpected UMI-tools version: {args.version!r}")
    if not args.sample or "\t" in args.sample or "\n" in args.sample:
        raise ValueError("sample identifier must be a nonempty single TSV field")

    policy = json.loads(args.policy.read_text(encoding="utf-8"))
    if (
        policy.get("schema_version") != 1
        or policy.get("representation") != "rx_v1"
        or policy.get("umi_length") != 8
        or policy.get("mapping_eligibility") != "all_primary_proper"
    ):
        raise ValueError("UMI preparation policy does not match the production contract")

    log = args.log.read_text(encoding="utf-8")
    if not log.strip():
        raise ValueError("UMI-tools log is empty")

    input_values = exactly_one(PATTERNS["tool_input_templates"], log, "input")
    if input_values[0] != input_values[1]:
        raise ValueError(f"UMI-tools input/read-pair counts disagree: {input_values}")
    parsed: dict[str, int | float] = {
        "tool_input_templates": int(input_values[0])
    }
    for name in (
        "molecules",
        "positions",
        "max_unique_umis_per_position",
        "unmatched_alignments",
        "mates_never_found",
    ):
        parsed[name] = int(exactly_one(PATTERNS[name], log, name)[0])
    parsed["mean_unique_umis_per_position"] = float(
        exactly_one(
            PATTERNS["mean_unique_umis_per_position"],
            log,
            "mean unique UMIs per position",
        )[0]
    )
    if not math.isfinite(parsed["mean_unique_umis_per_position"]):
        raise ValueError("mean UMI count is nonfinite")

    eligible = policy.get("eligible_acgt_templates")
    if isinstance(eligible, bool) or not isinstance(eligible, int) or eligible <= 0:
        raise ValueError("eligible UMI template denominator is invalid")
    if parsed["tool_input_templates"] != eligible:
        raise ValueError(
            "UMI-tools input differs from eligible policy denominator: "
            f"{parsed['tool_input_templates']} != {eligible}"
        )
    if parsed["mates_never_found"] != 0:
        raise ValueError(f"UMI-tools could not find {parsed['mates_never_found']} mates")
    molecules = parsed["molecules"]
    if isinstance(molecules, bool) or not isinstance(molecules, int):
        raise ValueError("UMI molecule count is invalid")
    if not 0 <= molecules <= eligible:
        raise ValueError("UMI molecule count is outside the eligible denominator")

    duplicates = eligible - molecules
    percentage = 100.0 * duplicates / eligible
    metrics = {
        "schema_version": 1,
        "algorithm": "umi_tools_directional",
        "acceptance_class": "scientific-truth",
        "representation": "rx_v1",
        "tool_version": args.version.strip(),
        "container": args.container,
        "umi_length": 8,
        "edit_distance": 1,
        "random_seed": 12345,
        "mapping_eligibility": "all_primary_proper",
        "paired_policy": {
            "unpaired": "discard",
            "chimeric": "discard",
            "unmapped": "discard",
            "multimapping_detection": "NH",
        },
        "umi_input_fragments": policy["primary_proper_templates"],
        "umi_primary_complete_templates": policy["primary_proper_templates"],
        "umi_acgt_primary_proper_templates": eligible,
        "umi_eligible_acgt_templates": eligible,
        "umi_invalid_n_templates": policy["excluded_n_umi_templates"],
        "umi_nh1_templates": policy["eligible_nh1_templates"],
        "umi_nh_gt1_templates": policy["eligible_multimapped_templates"],
        "umi_mapping_excluded_templates": 0,
        "umi_molecules": molecules,
        "umi_duplicate_templates": duplicates,
        "pct_umi_dup_eligible": percentage,
        **parsed,
    }
    args.metrics.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.report.write_text(
        "Sample\tpct_umi_dup\n"
        f"{args.sample}\t{rendered_percentage(percentage)}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
