#!/usr/bin/env python3
"""Write wdl/merge_cohort.wdl inputs from the per-sample outputs of one or more finished runs.

Output paths come from workflow metadata, recovered through Caper if the GCS export is missing,
so call-cache hits resolve to the earlier run directory that holds their files. Samples without
complete per-sample results are merged as failed: they are absent from the matrices and marked
status=failed in the sample sheet. Later run roots take precedence.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

from make_json_rnaseq import IMAGE_ROLES, MERGE_RESULTS_DOCKER, REFERENCE_ROLES

OUTPUTS = {
    task: {field: field + "_" + mode for field in ("genes", "isoforms", "log", "convergence", "gene_convergence")}
    for task, mode in (("umi_molecule_rsem", "umi"), ("rsem_quant", "all"))
}
OUTPUTS.update({
    "umi_molecule_feature_counts_task": {"fc_out": "fc_umi"}, "feature_counts": {"fc_out": "fc_all"},
    "qc_report": {"rnaseq_report": "qc", "diagnostics": "diagnostics"},
    "star_align": {"bam_file": "genomic_bam", "transcriptome_bam": "transcriptome_bam"},
    "udup": {"molecule_transcriptome_bam": "molecule_bam", "umi_metrics": "umi_metrics",
             "molecule_expression_metrics": "umi_expression_metrics"},
    "combined_contamination_qc": {"sampling_manifest": "sampling_manifest"},
    "trim_i1": {"trimmed_index": "trimmed_i1"},
})
HEADER = ["sample", "reference_release", "umi_available", "not_deduplicated", "expression_mode", "umi_status"]
BAMS = ("genomic_bam", "transcriptome_bam", "molecule_bam")
SCIENTIFIC_INPUTS = (REFERENCE_ROLES | IMAGE_ROLES | {"minimumLength", "index_adapter", "univ_adapter"}) - {
    "merge_results_docker", "multiqc_docker", "attach_umi_docker", "umi_dup_docker"}
DEFAULTS = {"reference_release": "unspecified", "run_pretrim_fastqc": True, "run_posttrim_fastqc": True,
            "run_contamination_qc": True, "combine_contamination_qc": False, "contamination_qc_pairs": 0,
            "run_alignment_qc": True, "run_umi_qc": True, "trim_trailing_i1_base": False,
            "use_umi_molecule_expression": True}


def workflow_metadata(root):
    root = root.rstrip("/")
    uri = root + "/metadata.json"
    result = subprocess.run(["gcloud", "storage", "cat", uri], capture_output=True, text=True)
    if result.returncode == 0:
        return json.loads(result.stdout)
    if "matched no objects" not in result.stderr and not re.search(r"\b404\b", result.stderr):
        raise RuntimeError(f"cannot read {uri}: {result.stderr.strip()}")

    workflow_id = root.rsplit("/", 1)[-1]
    print(f"Missing {uri}; retrieving metadata from Caper for {workflow_id}", file=sys.stderr)
    try:
        result = subprocess.run(["caper", "metadata", workflow_id], capture_output=True, text=True)
    except FileNotFoundError:
        raise RuntimeError("metadata recovery requires caper; run on the controller with access to this workflow") from None
    if result.returncode:
        raise RuntimeError(f"cannot recover metadata for {workflow_id}: {result.stderr.strip()}")
    metadata = json.loads(result.stdout)
    if metadata.get("id") != workflow_id or metadata.get("workflowRoot", "").rstrip("/") != root:
        raise ValueError(f"recovered metadata does not match run root {root}")
    return metadata


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


def resolved_inputs(inputs):
    prefix = "rnaseq_pipeline."
    return {**DEFAULTS, **{key[len(prefix):] if key.startswith(prefix) else key: value for key, value in inputs.items()}}


def sample_inputs(inputs):
    samples = inputs["sample_prefix"]
    indexes = inputs.get("fastq_index") or [""] * len(samples)
    arrays = (inputs["fastq1"], inputs["fastq2"], indexes)
    if not samples or len(set(samples)) != len(samples) or any(len(array) != len(samples) for array in arrays):
        raise ValueError("sample IDs must be unique and FASTQ arrays must match their length")
    return {sample: (inputs["fastq1"][i], inputs["fastq2"][i], indexes[i] or "") for i, sample in enumerate(samples)}


def check_compatible(requested, source, sample, raw, source_raw):
    if raw != source_raw:
        raise ValueError(f"raw FASTQs differ for {sample}")
    keys = SCIENTIFIC_INPUTS | (DEFAULTS.keys() - {"use_umi_molecule_expression", "trim_trailing_i1_base"})
    if raw[2]:
        keys |= {"attach_umi_docker", "umi_dup_docker", "trim_trailing_i1_base", "use_umi_molecule_expression"}
    missing = sorted(key for key in keys if key not in requested or key not in source)
    different = sorted(key for key in keys if requested.get(key) != source.get(key))
    if missing or different:
        raise ValueError(f"incompatible results for {sample}: missing={missing}, different={different}")


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
            for field, kind in OUTPUTS.get(task, {}).items():
                value = call.get("outputs", {}).get(field)
                if isinstance(value, list):
                    if len(value) > 1:
                        raise ValueError(f"expected at most one {field} for {sample}")
                    value = value[0] if value else None
                if value:
                    results.setdefault(sample, {})[kind] = value
    return results


def build(inputs, roots, read_metadata=workflow_metadata, list_sizes=gcloud_sizes):
    inputs = resolved_inputs(inputs)
    get = inputs.get
    samples = get("sample_prefix")
    raw = sample_inputs(inputs)
    raw_paths = {sample: [path for path in paths if path] for sample, paths in raw.items()}
    release = get("reference_release", "unspecified")
    trim = get("trim_trailing_i1_base", False)
    found = {}
    for root in roots:
        metadata = read_metadata(root)
        if metadata.get("status") not in ("Succeeded", "Failed", "Aborted"):
            raise ValueError(f"{root} must be a finished workflow")
        source = resolved_inputs(metadata["inputs"])
        source_raw = sample_inputs(source)
        metadata = {**metadata, "inputs": source}
        for sample, record in run_results(metadata).items():
            if sample not in raw:
                continue
            check_compatible(inputs, source, sample, raw[sample], source_raw[sample])
            mode = "umi" if get("use_umi_molecule_expression") and raw[sample][2] else "all"
            if all(k in record for k in (f"genes_{mode}", f"isoforms_{mode}", f"fc_{mode}", "qc", "diagnostics")):
                found[sample] = (mode, record)

    if not found:
        raise ValueError("no completed samples found under the given run roots")
    sizes = list_sizes(*[uri for sample, (_, record) in found.items()
                         for uri in raw_paths[sample] + [record[kind] for kind in BAMS if kind in record]])
    rows, failed, size_rows = [HEADER + (["i1_layout"] if trim else [])], [], []
    files = {key: [] for key in ("genes", "isoforms", "fc", "qc", "diagnostics")}
    supporting_files, source_rows = [], [["sample", "kind", "source_uri"]]
    for sample in samples:
        has_index = bool(raw[sample][2])
        mode, record = found.get(sample, (None, {}))
        if mode is None:
            failed.append(sample)
            umi = get("use_umi_molecule_expression", True) and has_index
        else:
            umi = mode == "umi"
            source_rows.extend([sample, kind, uri]
                               for kind, uri in zip(("raw_r1", "raw_r2", "raw_i1"), raw[sample]) if uri)
            for key, kind in (("genes", f"genes_{mode}"), ("isoforms", f"isoforms_{mode}"), ("fc", f"fc_{mode}"),
                              ("qc", "qc"), ("diagnostics", "diagnostics")):
                files[key].append(record[kind])
                source_rows.append([sample, key, record[kind]])
            for kind in (f"log_{mode}", f"convergence_{mode}", f"gene_convergence_{mode}",
                         "umi_metrics", "umi_expression_metrics", "sampling_manifest"):
                if kind in record:
                    supporting_files.append(record[kind])
                    source_rows.append([sample, kind, record[kind]])
            bam = lambda kind: sizes.get(record.get(kind))
            size_rows.append([sample, gib(*[sizes.get(path) for path in raw_paths[sample]]),
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
            prefix + "supporting_files": supporting_files, prefix + "source_rows": source_rows,
            prefix + "output_report_name": get("output_report_name"),
            prefix + "merge_results_docker": MERGE_RESULTS_DOCKER}


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
