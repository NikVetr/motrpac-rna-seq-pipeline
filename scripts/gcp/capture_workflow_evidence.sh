#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "ERROR: $*" >&2
    exit 1
}

if [[ "$#" -ne 3 ]]; then
    die "usage: $0 METADATA_JSON OUTPUT_DIRECTORY EXPECTED_REVISION"
fi

metadata="$1"
output_dir="$2"
expected_revision="$3"
[[ "$expected_revision" =~ ^[0-9a-f]{40}$ ]] ||
    die "expected revision must be a full lowercase Git commit"

for command_name in cp dirname find gcloud git jq mkdir mktemp mv rm sed sha256sum sort wc xargs
do
    command -v "$command_name" >/dev/null || die "missing command: $command_name"
done

[[ -f "$metadata" ]] || die "metadata file does not exist: $metadata"
[[ ! -e "$output_dir" ]] || die "output already exists: $output_dir"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)" ||
    die "capture script is not inside a Git working tree"
repo_revision="$(git -C "$repo_root" rev-parse --verify HEAD)" ||
    die "cannot resolve repository revision"
[[ "$repo_revision" == "$expected_revision" ]] ||
    die "capture revision does not match the submitted revision"
untracked_files="$(git -C "$repo_root" ls-files --others --exclude-standard)" ||
    die "cannot inspect untracked repository files"
repo_clean=true
if ! git -C "$repo_root" diff --quiet --ignore-submodules -- ||
    ! git -C "$repo_root" diff --cached --quiet --ignore-submodules -- ||
    [[ -n "$untracked_files" ]]
then
    repo_clean=false
fi

workflow_id="$(jq -er '
    .id |
    select(type == "string") |
    select(test("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"))
' "$metadata")" || die "metadata lacks a valid workflow ID"
jq -e '.calls | type == "object"' "$metadata" >/dev/null ||
    die "metadata lacks a calls object"

output_parent="$(dirname -- "$output_dir")"
mkdir -p -- "$output_parent"
stage_dir="$(mktemp -d "${output_parent}/.workflow-evidence-${workflow_id}.XXXXXX")"
cleanup() {
    rm -rf -- "$stage_dir"
}
trap cleanup EXIT HUP INT TERM

mkdir -p -- "${stage_dir}/batch-jobs" "${stage_dir}/task-streams"
mkdir -p -- "${stage_dir}/top-level-outputs"
if jq -e '
    .submittedFiles.options |
    fromjson |
    [
      paths(scalars) as $path |
      ($path[-1] | tostring | ascii_downcase) |
      select(test("token|secret|password|credential|private_key|service_account_json"))
    ] |
    length > 0
' "$metadata" >/dev/null
then
    die "submitted workflow options contain credential-like fields; refusing capture"
fi
cp -p -- "$metadata" "${stage_dir}/metadata.json"
jq -e '.submittedFiles.inputs | fromjson' "$metadata" \
    >"${stage_dir}/submitted-inputs.json" ||
    die "metadata does not contain valid submitted inputs"
jq -e '.submittedFiles.options | fromjson' "$metadata" \
    >"${stage_dir}/submitted-options.json" ||
    die "metadata does not contain valid submitted options"
jq -n \
    --arg revision "$repo_revision" \
    --arg expected_submission_revision "$expected_revision" \
    --argjson clean "$repo_clean" '
    {
      schema_version: 1,
      revision: $revision,
      expected_submission_revision: $expected_submission_revision,
      clean: $clean
    }
' >"${stage_dir}/repository.json"

describe_object() {
    local uri="$1"
    local raw_object
    [[ "$uri" =~ ^gs://[^[:space:]]+$ ]] || return 1
    raw_object="$(mktemp "${stage_dir}/.gcs-object.XXXXXX")"
    if ! gcloud storage objects describe "$uri" --format=json >"$raw_object"
    then
        rm -f -- "$raw_object"
        return 1
    fi
    if ! jq -ec --arg uri "$uri" '
        {
          uri: $uri,
          bucket: (.bucket // ""),
          name: (.name // ""),
          generation: ((.generation // "") | tostring),
          metageneration: ((.metageneration // "") | tostring),
          size_bytes: ((.size // .size_bytes // "") | tostring),
          md5_base64: (.md5_hash // .md5Hash // null),
          crc32c_base64: (.crc32c_hash // .crc32c // null),
          content_type: (.content_type // .contentType // null),
          created_utc: (.creation_time // .timeCreated // null),
          updated_utc: (.update_time // .updated // null),
          storage_class: (.storage_class // .storageClass // null)
        } |
        select(.bucket != "") |
        select(.name != "") |
        select(.generation | test("^[1-9][0-9]*$")) |
        select(.size_bytes | test("^[0-9]+$"))
    ' "$raw_object"
    then
        rm -f -- "$raw_object"
        return 1
    fi
    rm -f -- "$raw_object"
}

gcs_uris="${stage_dir}/submitted-gcs-uris.txt"
jq -er '
    [.. | strings | select(startswith("gs://"))] |
    unique |
    if length == 0 then error("submitted inputs contain no GCS objects") else .[] end
' "${stage_dir}/submitted-inputs.json" >"$gcs_uris" ||
    die "could not enumerate submitted GCS input objects"

object_records="${stage_dir}/input-object-metadata.jsonl"
: >"$object_records"
while IFS= read -r uri
do
    describe_object "$uri" >>"$object_records" ||
        die "cannot capture immutable metadata for submitted GCS object: $uri"
done <"$gcs_uris"

jq -s 'sort_by(.uri)' "$object_records" >"${stage_dir}/input-objects.json"
input_object_count="$(jq -er 'length | select(. > 0)' "${stage_dir}/input-objects.json")" ||
    die "submitted GCS input-object manifest is empty"
rm -f -- "$gcs_uris" "$object_records"

top_level_outputs="${stage_dir}/top-level-output-uris.tsv"
printf '%s\n' $'output_name\tindex\turi' >"$top_level_outputs"
jq -er '
    .outputs |
    if type != "object" then error("metadata outputs are not an object") else . end |
    to_entries[] as $output |
    [$output.value | .. | strings | select(startswith("gs://"))] |
    to_entries[] |
    [$output.key, (.key | tostring), .value] |
    @tsv
' "$metadata" >>"$top_level_outputs" ||
    die "could not enumerate top-level workflow output objects"

output_records="${stage_dir}/output-object-metadata.jsonl"
output_files="${stage_dir}/top-level-output-files.tsv"
: >"$output_records"
printf '%s\n' $'output_name\tindex\turi\tlocal_file' >"$output_files"
max_output_bytes=$((256 * 1024 * 1024))
while IFS=$'\t' read -r output_name output_index uri
do
    [[ "$output_name" =~ ^[A-Za-z0-9_.-]+$ ]] ||
        die "unsafe top-level output name: $output_name"
    [[ "$output_index" =~ ^[0-9]+$ ]] ||
        die "invalid top-level output index for $output_name"
    short_name="${output_name##*.}"
    case "$short_name" in
        rsem_genes_count|rsem_genes_tpm|rsem_genes_fpkm|feature_counts_file)
            ;;
        qc_report_file|contamination_sampling_manifests|umi_metrics)
            ;;
        multiqc_prealign_reports|multiqc_postalign_reports)
            ;;
        umi_molecule_expression_metrics|umi_molecule_rsem_genes_count)
            ;;
        umi_molecule_rsem_genes_tpm|umi_molecule_rsem_genes_fpkm|umi_molecule_feature_counts)
            ;;
        *)
            die "top-level output is outside the bounded capture allowlist: $output_name"
            ;;
    esac

    output_record="$(mktemp "${stage_dir}/.output-object.XXXXXX")"
    describe_object "$uri" >"$output_record" ||
        die "cannot capture immutable metadata for workflow output: $uri"
    output_size="$(jq -er '.size_bytes | tonumber' "$output_record")" ||
        die "workflow output lacks a valid size: $uri"
    (( output_size <= max_output_bytes )) ||
        die "workflow output exceeds the 256-MiB capture limit: $uri"
    jq -ec \
        --arg output_name "$output_name" \
        --argjson output_index "$output_index" \
        '. + {output_name: $output_name, output_index: $output_index}' \
        "$output_record" >>"$output_records"
    rm -f -- "$output_record"

    local_file="${short_name}.${output_index}.artifact"
    gcloud storage cp "$uri" "${stage_dir}/top-level-outputs/${local_file}" \
        >/dev/null || die "cannot copy bounded workflow output: $uri"
    printf '%s\t%s\t%s\t%s\n' \
        "$output_name" "$output_index" "$uri" "top-level-outputs/${local_file}" \
        >>"$output_files"
done < <(sed '1d' "$top_level_outputs")

jq -s 'sort_by(.output_name, .output_index)' "$output_records" \
    >"${stage_dir}/output-objects.json"
output_object_count="$(jq -er 'length | select(. > 0)' \
    "${stage_dir}/output-objects.json")" ||
    die "top-level workflow output manifest is empty"
rm -f -- "$output_records"

attempts="${stage_dir}/attempts.tsv"
printf '%s\n' $'call\tshard\tattempt\texecution_status\tpreemptible\tjob_id\tmonitoring\tstdout\tstderr' \
    >"$attempts"
jq -er '
    .calls | to_entries[] as $call |
    $call.value[] |
    [
      $call.key,
      (.shardIndex | tostring),
      (.attempt | tostring),
      (.executionStatus // ""),
      (.preemptible | tostring),
      .jobId,
      (.monitoringLog // ""),
      (.stdout // ""),
      (.stderr // "")
    ] | @tsv
' "$metadata" >>"$attempts" || die "could not enumerate call attempts"

missing="${stage_dir}/missing-artifacts.tsv"
printf '%s\n' $'kind\tcall\tshard\tattempt\turi_or_job' >"$missing"
missing_count=0

while IFS=$'\t' read -r call shard attempt execution_status preemptible job_id monitoring stdout stderr
do
    [[ "$call" =~ ^[A-Za-z0-9_.-]+$ ]] || die "unsafe call name in metadata: $call"
    [[ "$shard" =~ ^-?[0-9]+$ ]] || die "invalid shard for $call: $shard"
    [[ "$attempt" =~ ^[1-9][0-9]*$ ]] || die "invalid attempt for $call: $attempt"
    [[ "$job_id" =~ ^projects/([A-Za-z0-9_-]+)/locations/([A-Za-z0-9_-]+)/jobs/([A-Za-z0-9_-]+)$ ]] ||
        die "invalid Batch job resource for $call: $job_id"

    project="${BASH_REMATCH[1]}"
    location="${BASH_REMATCH[2]}"
    job_name="${BASH_REMATCH[3]}"
    call_name="${call#rnaseq_pipeline.}"
    stem="${call_name}.shard-${shard}.attempt-${attempt}"

    if ! gcloud batch jobs describe "$job_name" \
        --project="$project" \
        --location="$location" \
        --format=json >"${stage_dir}/batch-jobs/${stem}.json"
    then
        rm -f -- "${stage_dir}/batch-jobs/${stem}.json"
        printf 'batch_job\t%s\t%s\t%s\t%s\n' \
            "$call" "$shard" "$attempt" "$job_id" >>"$missing"
        missing_count=$((missing_count + 1))
    fi

    for stream_name in monitoring stdout stderr
    do
        case "$stream_name" in
            monitoring) uri="$monitoring" ;;
            stdout) uri="$stdout" ;;
            stderr) uri="$stderr" ;;
        esac
        if [[ -z "$uri" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "$stream_name" "$call" "$shard" "$attempt" "MISSING_URI" >>"$missing"
            missing_count=$((missing_count + 1))
        elif [[ "$uri" != gs://* ]]; then
            die "non-GCS $stream_name URI for $call attempt $attempt: $uri"
        elif ! gcloud storage cp "$uri" \
            "${stage_dir}/task-streams/${stem}.${stream_name}" >/dev/null
        then
            rm -f -- "${stage_dir}/task-streams/${stem}.${stream_name}"
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "$stream_name" "$call" "$shard" "$attempt" "$uri" >>"$missing"
            missing_count=$((missing_count + 1))
        fi
    done
done < <(sed '1d' "$attempts")

jq -n \
    --arg workflow_id "$workflow_id" \
    --arg workflow_status "$(jq -er '.status' "$metadata")" \
    --arg repository_revision "$repo_revision" \
    --argjson repository_clean "$repo_clean" \
    --argjson attempt_count "$(sed '1d' "$attempts" | wc -l)" \
    --argjson input_object_count "$input_object_count" \
    --argjson output_object_count "$output_object_count" \
    --argjson missing_artifact_count "$missing_count" '
    {
      schema_version: 1,
      workflow_id: $workflow_id,
      workflow_status: $workflow_status,
      repository_revision: $repository_revision,
      repository_clean: $repository_clean,
      attempt_count: $attempt_count,
      submitted_gcs_object_count: $input_object_count,
      top_level_output_object_count: $output_object_count,
      missing_artifact_count: $missing_artifact_count,
      complete: ($missing_artifact_count == 0)
    }
' >"${stage_dir}/capture-status.json"

(
    cd "$stage_dir"
    find . -type f ! -name evidence-manifest.sha256 -print0 |
        sort -z |
        xargs -0 sha256sum >evidence-manifest.sha256
)

mv -- "$stage_dir" "$output_dir"
stage_dir=""
trap - EXIT HUP INT TERM

if (( missing_count > 0 )); then
    die "evidence bundle is incomplete (${missing_count} missing artifacts): $output_dir"
fi
echo "Captured workflow evidence: $output_dir"
