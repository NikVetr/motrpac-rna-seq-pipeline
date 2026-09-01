# Usage example: python3 make_json_rnaseq.py -g gs://example/fastq_raw -o . -r batch7_qc_metrics -a human -v gencode_v39 -n 1
import argparse
import json
import os
import re
from pathlib import Path


R1_SUFFIX = "_R1.fastq.gz"
R2_SUFFIX = "_R2.fastq.gz"
I1_SUFFIX = "_I1.fastq.gz"
OUTPUT_REPORT_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
GENERATED_JSON_PATTERN = re.compile(r"set[1-9][0-9]*_rnaseq\.json")
IMMUTABLE_IMAGE_PATTERN = re.compile(r".+@sha256:[0-9a-f]{64}")
REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_ROLES = {
    "star_index",
    "gtf_file",
    "rsem_reference",
    "globin_genome_dir_tar",
    "rrna_genome_dir_tar",
    "phix_genome_dir_tar",
    "ref_flat",
}
IMAGE_ROLES = {
    "fastqc_docker",
    "attach_umi_docker",
    "cutadapt_docker",
    "star_docker",
    "feature_counts_docker",
    "rsem_docker",
    "bowtie_docker",
    "picard_docker",
    "umi_dup_docker",
    "samtools_docker",
    "collect_qc_docker",
    "merge_results_docker",
}
COMPATIBILITY_ROLES = {
    "star_index_builder",
    "star_runtime",
    "rsem_reference_builder",
    "rsem_runtime",
}
DEFAULT_RELEASE_MANIFESTS = {
    ("human", "gencode_v47"): REPO_ROOT
    / "config"
    / "release-profiles"
    / "human-gencode-v47.json",
}
SUPPORTED_REFERENCES = {
    ("rat", "rn6"),
    ("rat", "rn7"),
    ("rat", "rn8"),
    ("human", "gencode_v39"),
    ("human", "gencode_v47"),
}


def output_report_stem(value):
    if not isinstance(value, str):
        raise ValueError("output_report_name must be a filename-safe string")
    stem = re.sub(r"(?:\.csv)+$", "", value)
    if not OUTPUT_REPORT_PATTERN.fullmatch(stem) or stem in (".", ".."):
        raise ValueError("output_report_name must be a filename-safe stem or CSV name")
    return stem


def as_gcs_uri(path):
    return path if path.startswith("gs://") else "gs://" + path


def split_batches(values, count):
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError("num_chunks must be a positive integer")
    if not values:
        raise ValueError("no R1 FASTQs remain after filtering")
    if count > len(values):
        raise ValueError("num_chunks cannot exceed the filtered R1 FASTQ count")
    quotient, remainder = divmod(len(values), count)
    batches = []
    start = 0
    for index in range(count):
        size = quotient + (1 if index < remainder else 0)
        batches.append(values[start : start + size])
        start += size
    return batches


def build_batches(r1_paths, num_chunks, include_undetermined=False, include_index=False):
    selected = sorted(
        as_gcs_uri(path)
        for path in r1_paths
        if include_undetermined or "Undetermined_" not in os.path.basename(path)
    )
    if len(selected) != len(set(selected)):
        raise ValueError("R1 FASTQ listing contains duplicate objects")

    result = []
    seen_samples = set()
    for r1_batch in split_batches(selected, num_chunks):
        sample_prefix = []
        r2_batch = []
        i1_batch = []
        for r1 in r1_batch:
            basename = os.path.basename(r1)
            if not basename.endswith(R1_SUFFIX) or basename == R1_SUFFIX:
                raise ValueError("R1 object has an invalid filename: {}".format(r1))
            sample = basename[: -len(R1_SUFFIX)]
            if not OUTPUT_REPORT_PATTERN.fullmatch(sample) or sample in (".", ".."):
                raise ValueError("sample prefix is not filename-safe: {}".format(sample))
            if sample in seen_samples:
                raise ValueError("duplicate sample prefix: {}".format(sample))
            seen_samples.add(sample)
            sample_prefix.append(sample)
            r2_batch.append(r1[: -len(R1_SUFFIX)] + R2_SUFFIX)
            if include_index:
                i1_batch.append(r1[: -len(R1_SUFFIX)] + I1_SUFFIX)
        result.append(
            {
                "r1": r1_batch,
                "r2": r2_batch,
                "i1": i1_batch if include_index else None,
                "sample_prefix": sample_prefix,
            }
        )
    return result


def _require_exact_roles(profile, section_name, expected_roles):
    section = profile.get(section_name)
    if not isinstance(section, dict):
        raise ValueError("release manifest {} must be an object".format(section_name))
    actual_roles = set(section)
    missing = sorted(expected_roles - actual_roles)
    unknown = sorted(actual_roles - expected_roles)
    if missing or unknown:
        details = []
        if missing:
            details.append("missing {}".format(", ".join(missing)))
        if unknown:
            details.append("unknown {}".format(", ".join(unknown)))
        raise ValueError(
            "release manifest {} keys are invalid: {}".format(
                section_name, "; ".join(details)
            )
        )
    return section


def load_release_manifest(path, organism, version):
    manifest_path = Path(path)
    try:
        with manifest_path.open(encoding="utf-8") as file:
            profile = json.load(file)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(
            "cannot read release manifest {}: {}".format(manifest_path, exc)
        ) from exc

    if not isinstance(profile, dict):
        raise ValueError("release manifest root must be an object")
    if type(profile.get("schema_version")) is not int or profile["schema_version"] != 1:
        raise ValueError("release manifest schema_version must equal 1")
    profile_id = profile.get("profile_id")
    if not isinstance(profile_id, str) or not profile_id.strip():
        raise ValueError("release manifest profile_id must be a nonempty string")
    if profile.get("organism") != organism or profile.get("version") != version:
        raise ValueError(
            "release manifest {} declares {} {}, expected {} {}".format(
                profile_id,
                profile.get("organism"),
                profile.get("version"),
                organism,
                version,
            )
        )
    publication_state = profile.get("publication_state")
    if not isinstance(publication_state, str) or not publication_state.strip():
        raise ValueError(
            "release manifest publication_state must be a nonempty string"
        )

    references = _require_exact_roles(profile, "references", REFERENCE_ROLES)
    images = _require_exact_roles(profile, "images", IMAGE_ROLES)
    compatibility = _require_exact_roles(
        profile, "compatibility", COMPATIBILITY_ROLES
    )
    incomplete = sorted(
        "{}.{}".format(section_name, role)
        for section_name, section in (("references", references), ("images", images))
        for role, value in section.items()
        if not isinstance(value, str) or not value.strip() or value != value.strip()
    )
    if incomplete:
        if publication_state == "publication_pending":
            raise ValueError(
                "release profile {} is publication-pending; populate: {}".format(
                    profile_id, ", ".join(incomplete)
                )
            )
        raise ValueError(
            "release profile {} has empty execution values: {}".format(
                profile_id, ", ".join(incomplete)
            )
        )

    invalid_versions = sorted(
        role
        for role, value in compatibility.items()
        if not isinstance(value, str) or not value.strip() or value != value.strip()
    )
    if invalid_versions:
        raise ValueError(
            "release manifest compatibility values must be nonempty: {}".format(
                ", ".join(invalid_versions)
            )
        )
    if compatibility["star_index_builder"] != compatibility["star_runtime"]:
        raise ValueError(
            "release manifest STAR index-builder/runtime versions are incompatible: "
            "{} != {}".format(
                compatibility["star_index_builder"], compatibility["star_runtime"]
            )
        )
    if compatibility["rsem_reference_builder"] != compatibility["rsem_runtime"]:
        raise ValueError(
            "release manifest RSEM reference-builder/runtime versions are "
            "incompatible: {} != {}".format(
                compatibility["rsem_reference_builder"],
                compatibility["rsem_runtime"],
            )
        )

    mutable_images = sorted(
        role
        for role, value in images.items()
        if not IMMUTABLE_IMAGE_PATTERN.fullmatch(value)
    )
    if mutable_images:
        raise ValueError(
            "modern release manifest images must use immutable sha256 digests: {}".format(
                ", ".join(mutable_images)
            )
        )

    return {
        "rnaseq_pipeline.{}".format(role): value
        for section in (references, images)
        for role, value in section.items()
    }


def resolve_release_inputs(organism, version, release_manifest=None):
    manifest_path = release_manifest
    if manifest_path is None:
        manifest_path = DEFAULT_RELEASE_MANIFESTS.get((organism, version))
    if manifest_path is None:
        return None
    return load_release_manifest(manifest_path, organism, version)


def main(command_args: argparse.Namespace):
    if (command_args.organism, command_args.version) not in SUPPORTED_REFERENCES:
        raise ValueError(
            "unsupported organism/reference combination: {} {}".format(
                command_args.organism, command_args.version
            )
        )
    output_report_stem(command_args.output_report_name)
    use_umi_molecule_expression = getattr(
        command_args, "umi_molecule_expression", False
    )
    qc_settings = {
        "run_pretrim_fastqc": not getattr(command_args, "skip_pretrim_fastqc", False),
        "run_posttrim_fastqc": not getattr(command_args, "skip_posttrim_fastqc", False),
        "run_contamination_qc": not getattr(command_args, "skip_contamination_qc", False),
        "combine_contamination_qc": getattr(
            command_args, "combine_contamination_qc", False
        ),
        "contamination_qc_pairs": getattr(command_args, "contamination_qc_pairs", 0),
        "run_alignment_qc": not getattr(command_args, "skip_alignment_qc", False),
        "run_umi_qc": not getattr(command_args, "skip_umi_qc", False),
    }
    if use_umi_molecule_expression and not command_args.index:
        raise ValueError(
            "--umi-deduplicated-expression requires --index and matched I1 FASTQs"
        )
    release_inputs = resolve_release_inputs(
        command_args.organism,
        command_args.version,
        getattr(command_args, "release_manifest", None),
    )

    try:
        import gcsfs
    except ImportError as exc:
        raise RuntimeError("gcsfs is required to list the configured GCS input path") from exc

    if not command_args.gcp_path or not command_args.gcp_path.startswith("gs://"):
        raise ValueError("gcp_path must be a gs:// URI")
    output_path = Path(command_args.output_path)
    if not output_path.is_dir():
        raise ValueError("output_path must be an existing directory")
    existing_outputs = sorted(
        path.name
        for path in output_path.iterdir()
        if path.is_file() and GENERATED_JSON_PATTERN.fullmatch(path.name)
    )
    if existing_outputs:
        raise ValueError(
            "output_path already contains generated set JSONs: {}".format(
                ", ".join(existing_outputs)
            )
        )
    docker_repo = command_args.docker_repo.rstrip("/").strip()
    if not docker_repo:
        raise ValueError("docker_repo must be nonempty")

    fs = gcsfs.GCSFileSystem(project=command_args.project)
    batches = build_batches(
        fs.glob(command_args.gcp_path.rstrip("/") + "/*_R1.fastq.gz"),
        command_args.num_chunks,
        include_undetermined=command_args.undetermined,
        include_index=command_args.index,
    )
    expected_inputs = [
        path
        for batch in batches
        for key in ("r2", "i1")
        for path in (batch[key] or [])
    ]
    missing = [path for path in expected_inputs if not fs.exists(path)]
    if missing:
        raise ValueError(
            "required mate/index objects are missing: {}".format(", ".join(missing))
        )

    print("Number of batches to split:\t{}".format(command_args.num_chunks))
    documents = []
    for batch in batches:
        documents.append(
            make_json_dict(
                command_args.organism,
                command_args.version,
                docker_repo,
                command_args.output_report_name,
                batch["r1"],
                batch["r2"],
                batch["i1"],
                batch["sample_prefix"],
                release_inputs=release_inputs,
                use_umi_molecule_expression=use_umi_molecule_expression,
                **qc_settings,
            )
        )

    for batch_num, document in enumerate(documents, start=1):
        destination = output_path / "set{}_rnaseq.json".format(batch_num)
        temporary = destination.with_name(destination.name + ".tmp")
        try:
            with temporary.open("x", encoding="utf-8") as file:
                json.dump(obj=document, fp=file, indent=4)
                file.write("\n")
            os.replace(temporary, destination)
        finally:
            if temporary.exists():
                temporary.unlink()

    print("Success! Finished generating input jsons")
    return 0


def make_json_dict(
    organism,
    version,
    docker_repo,
    output_report_name,
    r1=None,
    r2=None,
    i1=None,
    prefix_list=None,
    release_inputs=None,
    use_umi_molecule_expression=False,
    run_pretrim_fastqc=True,
    run_posttrim_fastqc=True,
    run_contamination_qc=True,
    combine_contamination_qc=False,
    contamination_qc_pairs=0,
    run_alignment_qc=True,
    run_umi_qc=True,
):
    if r1 is None:
        r1 = []
    if r2 is None:
        r2 = []
    if prefix_list is None:
        prefix_list = []
    if not r1 or len(r1) != len(r2) or len(r1) != len(prefix_list):
        raise ValueError("R1, R2, and sample-prefix arrays must be nonempty and aligned")
    if i1 is not None and len(i1) != len(r1):
        raise ValueError("I1 array must be absent or aligned with R1 and R2")
    if use_umi_molecule_expression and not i1:
        raise ValueError(
            "UMI molecule expression requires a matched I1 FASTQ for every sample"
        )
    if type(contamination_qc_pairs) is not int or contamination_qc_pairs < 0:
        raise ValueError("contamination_qc_pairs must be a nonnegative integer")
    if not run_contamination_qc and (
        combine_contamination_qc or contamination_qc_pairs != 0
    ):
        raise ValueError(
            "combined or sampled contamination QC requires run_contamination_qc"
        )
    if len(prefix_list) != len(set(prefix_list)) or any(not value for value in prefix_list):
        raise ValueError("sample prefixes must be nonempty and unique")
    fastq_uris = r1 + r2 + (i1 or [])
    if len(fastq_uris) != len(set(fastq_uris)):
        raise ValueError("FASTQ URIs must be unique across R1, R2, and I1 roles")
    output_report_name = output_report_stem(output_report_name)

    if release_inputs is not None:
        organism_references = {}
    elif organism == "rat" and version == "rn6":
        organism_references = {
            "rnaseq_pipeline.star_index": "gs://omicspipelines-public-resources/rnaseq/references/rat/Rnor6_v96_star_index.tar.gz",
            "rnaseq_pipeline.gtf_file": "gs://omicspipelines-public-resources/rnaseq/references/rat/Rattus_norvegicus.Rnor_6.0.96.gtf",
            "rnaseq_pipeline.rsem_reference": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn6_rsem_reference.tar.gz",
            "rnaseq_pipeline.globin_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_globin.tar.gz",
            "rnaseq_pipeline.rrna_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_rRNA.tar.gz",
            "rnaseq_pipeline.phix_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/phix.tar.gz",
            "rnaseq_pipeline.ref_flat": "gs://omicspipelines-public-resources/rnaseq/references/rat/refFlat_rn6_v96.txt",
        }
    elif organism == "rat" and version == "rn7" :
        organism_references = {
            "rnaseq_pipeline.star_index": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/rn7_v108_star_index.tar.gz",
            "rnaseq_pipeline.gtf_file": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/Rattus_norvegicus.mRatBN7.2.108.gtf",
            "rnaseq_pipeline.rsem_reference": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/rn7_rsem_reference.tar.gz",
            "rnaseq_pipeline.globin_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_globin.tar.gz",
            "rnaseq_pipeline.rrna_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_rRNA.tar.gz",
            "rnaseq_pipeline.phix_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/phix.tar.gz",
            "rnaseq_pipeline.ref_flat": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/refFlat_mRatBN7.2_v108.txt",
        }
    elif organism == "rat" and version == "rn8":
        organism_references = {
            "rnaseq_pipeline.star_index": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/rn8_v115_star_index.tar.gz",
            "rnaseq_pipeline.gtf_file": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/Rattus_norvegicus.GRCr8.115.gtf",
            "rnaseq_pipeline.rsem_reference": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/rn8_rsem_reference.tar.gz",
            "rnaseq_pipeline.globin_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_globin.tar.gz",
            "rnaseq_pipeline.rrna_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn_rRNA.tar.gz",
            "rnaseq_pipeline.phix_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/rat/phix.tar.gz",
            "rnaseq_pipeline.ref_flat": "gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/refFlat_GRCr8_v115.txt",
        }
    elif organism == "human" and version == "gencode_v39" :
        organism_references = {
            "rnaseq_pipeline.star_index": "gs://omicspipelines-public-resources/rnaseq/references/human/hg38_v39_star_index.tar.gz",
            "rnaseq_pipeline.gtf_file": "gs://omicspipelines-public-resources/rnaseq/references/human/GRCh38.v39.primary_assembly.annotation.gtf",
            "rnaseq_pipeline.rsem_reference": "gs://omicspipelines-public-resources/rnaseq/references/human/hg38_rsem_reference.tar.gz",
            "rnaseq_pipeline.globin_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/hs_globin.tar.gz",
            "rnaseq_pipeline.rrna_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/hs_rRNA.tar.gz",
            "rnaseq_pipeline.phix_genome_dir_tar": "gs://omicspipelines-public-resources/rnaseq/references/human/phix.tar.gz",
            "rnaseq_pipeline.ref_flat": "gs://omicspipelines-public-resources/rnaseq/references/human/refFlat_hg38_v39.txt",
        }
    else:
        raise ValueError("unsupported organism/reference combination: {} {}".format(organism, version))

    filled_dict = {
        "rnaseq_pipeline.fastq1": r1,
        "rnaseq_pipeline.fastq2": r2,
        "rnaseq_pipeline.fastq_index": i1,
        "rnaseq_pipeline.sample_prefix": prefix_list,
        "rnaseq_pipeline.pretrim_fastqc_ncpu": 8,
        "rnaseq_pipeline.pretrim_fastqc_ramGB": 40,
        "rnaseq_pipeline.pretrim_fastqc_disk": 100,
        "rnaseq_pipeline.fastqc_docker": f"{docker_repo}/fastqc:latest",
        "rnaseq_pipeline.attach_umi_ncpu": 8,
        "rnaseq_pipeline.attach_umi_ramGB": 40,
        "rnaseq_pipeline.attach_umi_disk": 100,
        "rnaseq_pipeline.attach_umi_docker": f"{docker_repo}/umi_attach:latest",
        "rnaseq_pipeline.minimumLength": 20,
        "rnaseq_pipeline.index_adapter": "AGATCGGAAGAGC",
        "rnaseq_pipeline.univ_adapter": "AGATCGGAAGAGC",
        "rnaseq_pipeline.cutadapt_ncpu": 8,
        "rnaseq_pipeline.cutadapt_ramGB": 45,
        "rnaseq_pipeline.cutadapt_disk": 100,
        "rnaseq_pipeline.cutadapt_docker": f"{docker_repo}/cutadapt:latest",
        "rnaseq_pipeline.posttrim_fastqc_ncpu": 8,
        "rnaseq_pipeline.posttrim_fastqc_ramGB": 36,
        "rnaseq_pipeline.posttrim_fastqc_disk": 100,
        "rnaseq_pipeline.star_ncpu": 12,
        "rnaseq_pipeline.star_ramGB": 120,
        "rnaseq_pipeline.star_disk": 400,
        "rnaseq_pipeline.star_docker": f"{docker_repo}/star:latest",
        "rnaseq_pipeline.feature_counts_ncpu": 8,
        "rnaseq_pipeline.feature_counts_ramGB": 48,
        "rnaseq_pipeline.feature_counts_disk": 100,
        "rnaseq_pipeline.feature_counts_docker": f"{docker_repo}/feature_counts:latest",
        "rnaseq_pipeline.rsem_ncpu": 10,
        "rnaseq_pipeline.rsem_ramGB": 48,
        "rnaseq_pipeline.rsem_disk": 150,
        "rnaseq_pipeline.rsem_docker": f"{docker_repo}/rsem:latest",
        "rnaseq_pipeline.bowtie2_globin_ncpu": 12,
        "rnaseq_pipeline.bowtie2_globin_ramGB": 80,
        "rnaseq_pipeline.bowtie2_globin_disk": 200,
        "rnaseq_pipeline.bowtie2_rrna_ncpu": 12,
        "rnaseq_pipeline.bowtie2_rrna_ramGB": 80,
        "rnaseq_pipeline.bowtie2_rrna_disk": 200,
        "rnaseq_pipeline.bowtie2_phix_ncpu": 12,
        "rnaseq_pipeline.bowtie2_phix_ramGB": 80,
        "rnaseq_pipeline.bowtie2_phix_disk": 200,
        "rnaseq_pipeline.bowtie_docker": f"{docker_repo}/bowtie:latest",
        "rnaseq_pipeline.markdup_ncpu": 10,
        "rnaseq_pipeline.markdup_ramGB": 96,
        "rnaseq_pipeline.markdup_disk": 300,
        "rnaseq_pipeline.rnaqc_ncpu": 10,
        "rnaseq_pipeline.rnaqc_ramGB": 48,
        "rnaseq_pipeline.rnaqc_disk": 100,
        "rnaseq_pipeline.picard_docker": f"{docker_repo}/picard:latest",
        "rnaseq_pipeline.umi_dup_ncpu": 8,
        "rnaseq_pipeline.umi_dup_ramGB": 36,
        "rnaseq_pipeline.umi_dup_disk": 200,
        "rnaseq_pipeline.umi_dup_docker": f"{docker_repo}/umi_dup:latest",
        "rnaseq_pipeline.mapped_ncpu": 8,
        "rnaseq_pipeline.mapped_ramGB": 36,
        "rnaseq_pipeline.mapped_disk": 200,
        "rnaseq_pipeline.samtools_docker": f"{docker_repo}/samtools:latest",
        "rnaseq_pipeline.collect_qc_ncpu": 8,
        "rnaseq_pipeline.collect_qc_ramGB": 16,
        "rnaseq_pipeline.collect_qc_disk": 100,
        "rnaseq_pipeline.collect_qc_docker": f"{docker_repo}/collect_qc:latest",
        "rnaseq_pipeline.output_report_name": output_report_name,
        "rnaseq_pipeline.merge_results_ncpu": 4,
        "rnaseq_pipeline.merge_results_ramGB": 16,
        "rnaseq_pipeline.merge_results_disk": 200,
        "rnaseq_pipeline.merge_results_docker": f"{docker_repo}/merge_results:latest",
    }
    if use_umi_molecule_expression:
        filled_dict["rnaseq_pipeline.use_umi_molecule_expression"] = True
    for name, enabled in {
        "run_pretrim_fastqc": run_pretrim_fastqc,
        "run_posttrim_fastqc": run_posttrim_fastqc,
        "run_contamination_qc": run_contamination_qc,
        "run_alignment_qc": run_alignment_qc,
        "run_umi_qc": run_umi_qc,
    }.items():
        if not enabled:
            filled_dict["rnaseq_pipeline.{}".format(name)] = False
    if combine_contamination_qc:
        filled_dict["rnaseq_pipeline.combine_contamination_qc"] = True
    if contamination_qc_pairs > 0:
        filled_dict["rnaseq_pipeline.contamination_qc_pairs"] = contamination_qc_pairs

    d = {**filled_dict, **organism_references, **(release_inputs or {})}

    return d


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script is used to generate input json files from the "
        "fastq_raw dir on gcp for running rna-seq pipeline on GCP "
    )
    parser.add_argument(
        "-g",
        "--gcp_path",
        help="location of the submission batch directory in gcp that contains the "
        "fastq_raw dir",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        help="output path, where you want the input jsons to be written",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r",
        "--output_report_name",
        help="suffix-free report stem or name ending in one or more .csv suffixes",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-u",
        "--undetermined",
        help="Adding this flag will process undetermined FastQ files if they exist. "
        'These are fastq files with prefix "Undetermined_". If this flag isn\'t '
        'passed, items with prefix "Undetermined_" will be removed',
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--organism",
        help="organism name, e.g. rat or human",
        choices=["rat", "human"],
        required=True,
    )
    parser.add_argument(
        "-v",
        "--version",
        help="genome build version to use for references",
        choices=["rn6", "rn7", "rn8", "gencode_v39", "gencode_v47"],
        required=True,
    )
    parser.add_argument(
        "-n",
        "--num_chunks",
        help="number of chunks to split the input files, should always be <= number of "
        "input files",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-d",
        "--docker_repo",
        help="Docker repository prefix containing the images used in the workflow",
        type=str,
        default="us-docker.pkg.dev/motrpac-portal/rnaseq",
    )
    parser.add_argument(
        "--release-manifest",
        help="complete release profile overriding reference and Docker workflow inputs",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--index",
        help="Adding this flag will add index files to the input JSON",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--umi-deduplicated-expression",
        dest="umi_molecule_expression",
        help="add directional UMI molecule RSEM and featureCounts matrices; requires -i",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--skip-pretrim-fastqc",
        help="skip raw-read FastQC; the QC matrix retains explicit empty fields",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--skip-posttrim-fastqc",
        help="skip trimmed-read FastQC and leave its QC fields empty",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--skip-contamination-qc",
        help="skip globin, rRNA, and PhiX Bowtie2 screens",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--combine-contamination-qc",
        help="run full-depth globin, rRNA, and PhiX screens serially on one worker",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--contamination-qc-pairs",
        help=(
            "deterministically sample this many post-trim pairs for the combined "
            "screens; 0 retains full depth"
        ),
        type=int,
        default=0,
    )
    parser.add_argument(
        "--skip-alignment-qc",
        help="skip Picard duplicate/RNA metrics and chromosome summaries",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--skip-umi-qc",
        help="skip directional UMI QC unless molecule expression requires grouping",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-p", 
	"--project", 
	help="Project name on the google cloud platform", 
	type=str
    )
    args = parser.parse_args()
    raise SystemExit(main(args))
