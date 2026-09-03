#!/usr/bin/env bash

set -uo pipefail

project="${1:-motrpac-portal}"
batch_location="${2:-us-west1}"
controller_name="${3:-omicspipelines-get}"
controller_zone="${4:-us-west1-a}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
failures=0

section() {
    printf '\n## %s\n' "$1"
}

run_check() {
    local label="$1"
    local status
    shift

    if "$@"
    then
        return 0
    else
        status=$?
        failures=$((failures + 1))
        printf 'CHECK FAILED [%s] exit=%d; continuing read-only audit\n' \
            "$label" "$status" >&2
    fi
}

controller_startup_summary() {
    gcloud compute instances describe "$controller_name" \
        --project="$project" \
        --zone="$controller_zone" \
        --format=json \
        | jq '[.metadata.items[]? | select(.key == "startup-script" or .key == "startup-script-url") | .value | test("cromwell|caper|docker|mysql|systemctl|supervisor"; "i")] | {startup_script_count: length, mentions_pipeline_services: any}'
}

controller_configuration() {
    gcloud compute instances describe "$controller_name" \
        --project="$project" \
        --zone="$controller_zone" \
        --format=json \
        | jq '{
            name,
            status,
            zone: (.zone | split("/")[-1]),
            machineType: (.machineType | split("/")[-1]),
            deletionProtection,
            scheduling,
            labels,
            serviceAccounts,
            metadata_keys: [.metadata.items[]?.key],
            disks: [.disks[]? | {
                deviceName,
                boot,
                autoDelete,
                source: (.source | split("/")[-1])
            }]
        }'
}

describe_bucket() {
    gcloud storage buckets describe "$1" --raw --format=json \
        | jq '{name,location,locationType,storageClass,iamConfiguration,retentionPolicy}'
}

regional_quota_summary() {
    gcloud compute regions describe "$1" \
        --project="$project" \
        --format=json \
        | jq -r '
            ["METRIC", "LIMIT", "USAGE"],
            (.quotas[]
                | select(.metric | test("^(CPUS|PREEMPTIBLE_CPUS|SPOT_CPUS|E2_CPUS|N2_CPUS|N2D_CPUS|C2_CPUS|SSD_TOTAL_GB|IN_USE_ADDRESSES)$"))
                | [.metric, .limit, .usage])
            | @tsv
        '
}

cd "$repo_root"

section "Modernization revision"
git status --short --branch
git rev-parse HEAD

section "Active gcloud identity and defaults"
run_check "gcloud version" gcloud version
run_check "active gcloud identity" \
    gcloud auth list --filter='status:ACTIVE' --format='table(account,status)'
run_check "gcloud defaults" gcloud config list \
    --format='yaml(core.account,core.project,compute.region,compute.zone,batch.location)'
run_check "project metadata" \
    gcloud projects describe "$project" --format='value(projectId,projectNumber)'
if gcloud auth application-default print-access-token > /dev/null 2>&1
then
    printf 'Application Default Credentials: available\n'
else
    printf 'Application Default Credentials: missing\n'
fi

section "Existing controller (read only; it remains stopped)"
run_check "controller configuration" controller_configuration
run_check "controller startup metadata" controller_startup_summary

section "Candidate pipeline service accounts"
run_check "pipeline service-account inventory" \
    gcloud iam service-accounts list \
    --project="$project" \
    --filter='email:batch OR email:cromwell OR email:pipeline OR email:compute' \
    --format='table(email,displayName,disabled)'

section "Active Batch jobs and running VMs"
run_check "active Batch jobs" \
    gcloud batch jobs list \
    --project="$project" \
    --location="$batch_location" \
    --filter='status.state=(QUEUED,SCHEDULED,RUNNING)' \
    --format='table(name,status.state,createTime)'
run_check "running Compute Engine VMs" \
    gcloud compute instances list \
    --project="$project" \
    --filter='status=RUNNING' \
    --format='table(name,zone.basename(),machineType.basename(),status,labels)'
run_check "Batch worker disks" \
    gcloud compute disks list \
    --project="$project" \
    --filter='labels.batch-node:*' \
    --format='table(name,zone.basename(),sizeGb,status,labels)'

section "Required enabled APIs"
run_check "enabled APIs" \
    gcloud services list \
    --enabled \
    --project="$project" \
    --filter='config.name=(batch.googleapis.com,compute.googleapis.com,storage.googleapis.com,artifactregistry.googleapis.com,cloudbuild.googleapis.com,logging.googleapis.com,monitoring.googleapis.com)' \
    --format='value(config.name)'

section "Existing Artifact Registry"
run_check "RNA-seq Artifact Registry" \
    gcloud artifacts repositories describe rnaseq \
    --project="$project" \
    --location=us \
    --format='yaml(name,format,mode,location,createTime)'

section "Relevant bucket locations and protections"
for bucket in \
    gs://omicspipelines-get \
    gs://motrpac-data-raw-human-main \
    gs://omicspipelines-public-resources
do
    run_check "bucket metadata: $bucket" describe_bucket "$bucket"
done

section "Representative input and reference access"
run_check "representative muscle FASTQ" \
    gcloud storage objects describe \
    gs://motrpac-data-raw-human-main/human-main/transcriptomics/t10-muscle/rna-seq/fastq/tr04/13146031002_R1.fastq.gz \
    --format='json(name,size,generation,md5_hash,crc32c,storage_class)'
run_check "representative RNA-seq reference" \
    gcloud storage objects describe \
    gs://omicspipelines-public-resources/rnaseq/references/human/hs_globin.tar.gz \
    --format='json(name,size,generation,md5_hash,crc32c,storage_class)'

section "Regional safety quotas"
for region in "$batch_location" us-west1
do
    run_check "regional quotas: $region" regional_quota_summary "$region"
done

printf '\nREAD_ONLY_PREFLIGHT_FAILURES=%d\n' "$failures"
printf 'Preflight completed without starting a VM, publishing an artifact, or submitting a job.\n'
if ((failures > 0))
then
    exit 1
fi
