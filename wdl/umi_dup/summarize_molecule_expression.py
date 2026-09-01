#!/usr/bin/env python3
"""Bind UMI grouping and propagated-alignment denominators into one record."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--umi-metrics", type=Path, required=True)
    parser.add_argument("--propagation-metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    umi = json.loads(args.umi_metrics.read_text(encoding="utf-8"))
    propagation = json.loads(
        args.propagation_metrics.read_text(encoding="utf-8")
    )
    if umi.get("algorithm") != "umi_tools_directional":
        raise ValueError("molecule expression requires directional UMI-tools grouping")
    if umi.get("acceptance_class") != "scientific-truth":
        raise ValueError("unexpected UMI grouping acceptance class")
    if propagation.get("status") != "production-shadow":
        raise ValueError("unexpected propagation status")
    qnames = propagation["qnames"]
    selected = qnames["selected_representative_qnames"]
    if selected != umi.get("umi_molecules"):
        raise ValueError(
            "representative QNAME and UMI molecule counts differ: "
            f"{selected} != {umi.get('umi_molecules')}"
        )
    genomic_records = propagation["genomic"]["kept_alignment_records"]
    transcriptome_records = propagation["transcriptome"][
        "kept_alignment_records"
    ]
    if genomic_records != qnames["selected_genomic_alignment_records"]:
        raise ValueError("featureCounts propagated-record denominator disagrees")
    if transcriptome_records != qnames[
        "selected_transcriptome_alignment_records"
    ]:
        raise ValueError("RSEM propagated-record denominator disagrees")

    result = {
        "schema_version": 1,
        "acceptance_class": "scientific-truth",
        "status": "production-shadow",
        "algorithm": "umi_tools_directional_molecule_expression_shadow_v1",
        "provenance": {
            "umi_algorithm": umi["algorithm"],
            "umi_tools_version": umi["tool_version"],
            "umi_container": umi["container"],
            "edit_distance": umi["edit_distance"],
            "random_seed": umi["random_seed"],
            "representation": umi["representation"],
            "propagation_algorithm": propagation["algorithm"],
            "propagation_runtime_versions": propagation["runtime_versions"],
        },
        "denominators": {
            "umi_eligible_genomic_templates": umi[
                "umi_eligible_acgt_templates"
            ],
            "selected_genomic_representative_qnames": selected,
            "selected_genomic_alignment_records_for_featurecounts": (
                genomic_records
            ),
            "transcriptome_source_unique_qnames": qnames[
                "transcriptome_source_unique_qnames"
            ],
            "selected_representative_qnames_present_in_transcriptome": qnames[
                "selected_transcriptome_present_qnames"
            ],
            "selected_representative_qnames_absent_from_transcriptome": qnames[
                "selected_transcriptome_absent_qnames"
            ],
            "selected_transcriptome_alignment_records_for_rsem": (
                transcriptome_records
            ),
        },
        "interpretation": {
            "conventional_outputs": "unchanged and remain authoritative",
            "featurecounts_shadow": (
                "one selected genomic UMI representative QNAME, with every "
                "original genomic alignment retained"
            ),
            "rsem_shadow": (
                "the transcript-compatible subset of those same representatives, "
                "with every STAR transcript alternative retained"
            ),
            "warning": (
                "featureCounts and RSEM molecule denominators differ when selected "
                "genomic representatives lack a compatible STAR transcript alignment"
            ),
        },
        "propagation": propagation,
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
