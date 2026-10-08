#!/usr/bin/env python3
"""Write wdl/merge_cohort.wdl inputs from the per-sample outputs of one or more finished runs.

Output paths come from the metadata.json that Caper writes to each workflow root, so call-cache
hits resolve to the earlier run directory that holds their files. Samples without complete
per-sample results are merged as failed: they are absent from the matrices and marked
status=failed in the sample sheet. Later run roots take precedence.
"""

import argparse
import json
import re
import subprocess
from pathlib import Path


KINDS = {  # (task, filename suffix after the sample ID) -> kind
    ("umi_molecule_rsem", ".genes.results"): "genes_umi", ("rsem_quant", ".genes.results"): "genes_all",
    ("umi_molecule_rsem", ".isoforms.results"): "isoforms_umi", ("rsem_quant", ".isoforms.results"): "isoforms_all",
    ("umi_molecule_feature_counts_task", ".out"): "fc_umi", ("feature_counts", ".out"): "fc_all",
    ("qc_report", "_qc_info.csv"): "qc", ("qc_report", ".qc_diagnostics.json"): "diagnostics",
    ("star_align", ".Aligned.sortedByCoord.out.bam"): "genomic_bam",
    ("star_align", ".Aligned.toTranscriptome.out.bam"): "transcriptome_bam",
    ("udup", ".umi_molecules.transcriptome.bam"): "molecule_bam", ("trim_i1", "_I1.trimmed.fastq.gz"): "trimmed_i1",
}
HEADER = ["sample", "reference_release", "umi_available", "not_deduplicated", "expression_mode", "umi_status"]
BAMS = ("genomic_bam", "transcriptome_bam", "molecule_bam")


def gcloud_metadata(root):
    return json.loads(subprocess.run(["gcloud", "storage", "cat", root.rstrip("/") + "/metadata.json"],
                                     check=True, capture_output=True, text=True).stdout)


def gcloud_sizes(*urls):
    """Return {uri: bytes} for the listed objects that exist; deleted intermediates are absent."""
    sizes = {}
    for start in range(0, len(urls), 500):
        listing = subprocess.run(["gcloud", "storage", "ls", "-l", *urls[start:start + 500]],
                                 capture_output=True, text=True)
        if listing.returncode and "matched no objects" not in listing.stderr:
            raise RuntimeError(listing.stderr)
        sizes.update({m[2]: int(m[1]) for m in re.finditer(r"^\s*(\d+)\s+\S+\s+(gs://\S+)$", listing.stdout, re.M)})
    return sizes


def gib(*sizes):
    return "" if None in sizes else "{:.6f}".format(sum(sizes) / 2**30)


def uris(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from uris(item)


def run_results(metadata):
    """Return {sample: {kind: uri}} from the completed calls of one workflow."""
    samples = metadata["inputs"]["sample_prefix"]
    results = {}
    for call_name, calls in metadata["calls"].items():
        task = call_name.rsplit(".", 1)[-1]
        for call in calls:
            if call["executionStatus"] != "Done" or call["shardIndex"] < 0:
                continue
            sample = samples[call["shardIndex"]]
            for uri in uris(list(call.get("outputs", {}).values())):
                filename = uri.rsplit("/", 1)[-1]
                kind = KINDS.get((task, filename[len(sample):])) if filename.startswith(sample) else None
                if kind:
                    results.setdefault(sample, {})[kind] = uri
    return results


def build(inputs, roots, read_metadata=gcloud_metadata, list_sizes=gcloud_sizes):
    get = lambda name, default=None: inputs.get("rnaseq_pipeline." + name, inputs.get(name, default))
    samples = get("sample_prefix")
    indexes = [path or "" for path in (get("fastq_index") or [""] * len(samples))]
    raw_paths = {sample: [path for path in (get("fastq1")[i], get("fastq2")[i], indexes[i]) if path]
                 for i, sample in enumerate(samples)}
    release = get("reference_release", "unspecified")
    trim = get("trim_trailing_i1_base", False)
    found = {}
    for root in roots:
        metadata = read_metadata(root)
        if metadata["inputs"].get("reference_release", "unspecified") != release:
            raise ValueError(f"{root} did not use reference_release {release}")
        for sample, record in run_results(metadata).items():
            mode = "umi" if "genes_umi" in record else "all"
            if sample in raw_paths and all(k in record for k in (f"genes_{mode}", f"isoforms_{mode}", f"fc_{mode}", "qc", "diagnostics")):
                found[sample] = (mode, record)

    if not found:
        raise ValueError("no completed samples found under the given run roots")
    sizes = list_sizes(*[uri for sample, (_, record) in found.items()
                         for uri in raw_paths[sample] + [record[kind] for kind in BAMS if kind in record]])
    missing = [path for sample in found for path in raw_paths[sample] if path not in sizes]
    if missing:
        raise ValueError(f"raw FASTQs not found: {missing}")
    rows, failed, size_rows = [HEADER + (["i1_layout"] if trim else [])], [], []
    files = {key: [] for key in ("genes", "isoforms", "fc", "qc", "diagnostics")}
    for shard, sample in enumerate(samples):
        has_index = bool(indexes[shard])
        mode, record = found.get(sample, (None, {}))
        if mode is None:
            failed.append(sample)
            umi = get("use_umi_molecule_expression", True) and has_index
        else:
            umi = mode == "umi"
            for key, kind in (("genes", f"genes_{mode}"), ("isoforms", f"isoforms_{mode}"), ("fc", f"fc_{mode}"),
                              ("qc", "qc"), ("diagnostics", "diagnostics")):
                files[key].append(record[kind])
            bam = lambda kind: sizes.get(record.get(kind))
            size_rows.append([sample, gib(*[sizes[path] for path in raw_paths[sample]]),
                              gib(bam("genomic_bam")), gib(bam("transcriptome_bam")),
                              gib(bam("molecule_bam") if umi else bam("transcriptome_bam"))])
        status = "skipped_no_umi" if not has_index else "deduplicated" if umi else "not_requested"
        row = [sample, release, "1" if has_index else "0", "0" if umi else "1",
               "umi_molecules" if umi else "all_read", status]
        if trim:
            row.append("none" if not has_index else "umi8_trailing_base_trimmed" if "trimmed_i1" in record else "umi8")
        rows.append(row)
    prefix = "rnaseq_merge."
    return {prefix + "sample_prefix": [sample for sample in samples if sample in found],
            prefix + "expression_metadata_rows": rows, prefix + "failed_samples": failed,
            prefix + "rsem_gene_results": files["genes"], prefix + "rsem_isoform_results": files["isoforms"],
            prefix + "feature_counts_files": files["fc"], prefix + "qc_report_files": files["qc"],
            prefix + "qc_diagnostics": files["diagnostics"], prefix + "sample_size_rows": size_rows,
            prefix + "output_report_name": get("output_report_name"),
            prefix + "merge_results_docker": get("merge_results_docker")}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="the cohort's submitted rnaseq_pipeline inputs JSON")
    parser.add_argument("--run-root", action="append", required=True,
                        help="gs://.../rnaseq_pipeline/WORKFLOW_ID; repeat for recovery runs, latest last")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    document = build(json.loads(args.inputs.read_text()), args.run_root)
    with args.output.open("x") as handle:
        json.dump(document, handle, indent=2)
        handle.write("\n")
    failed = document["rnaseq_merge.failed_samples"]
    print(f"{len(document['rnaseq_merge.sample_prefix'])} completed samples; {len(failed)} failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
