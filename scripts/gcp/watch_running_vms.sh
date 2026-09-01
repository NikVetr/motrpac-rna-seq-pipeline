#!/usr/bin/env bash

set -euo pipefail

project="${1:-motrpac-portal}"
max_new_vms="${2:-5}"
interval_seconds="${3:-15}"

if [[ ! "$max_new_vms" =~ ^[0-9]+$ ]]
then
    printf 'max_new_vms must be a nonnegative integer: %s\n' "$max_new_vms" >&2
    exit 2
fi
if [[ ! "$interval_seconds" =~ ^[1-9][0-9]*$ ]]
then
    printf 'interval_seconds must be a positive integer: %s\n' \
        "$interval_seconds" >&2
    exit 2
fi

list_running_vms() {
    gcloud compute instances list \
        --project="$project" \
        --filter='status=RUNNING' \
        --format=json
}

baseline_json="$(list_running_vms)"
baseline_names="$(jq '[.[].name]' <<< "$baseline_json")"
baseline_count="$(jq 'length' <<< "$baseline_json")"

printf 'Baseline captured: project=%s running_vms=%s max_new_vms=%s interval=%ss\n' \
    "$project" "$baseline_count" "$max_new_vms" "$interval_seconds"
printf 'Press Ctrl-C to stop. This monitor does not change cloud state.\n'

while true
do
    if ! current_json="$(list_running_vms)"
    then
        printf 'VM inventory failed at %s; retrying.\n' \
            "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
        sleep "$interval_seconds"
        continue
    fi

    counts="$(
        jq -r --argjson baseline "$baseline_names" '
            [
                length,
                ([.[] | select(.name as $name | ($baseline | index($name)) == null)] | length),
                ([.[]
                    | select(
                        (.labels // {}) as $labels
                        | ($labels | has("batch-node"))
                            or ($labels | has("batch-job-id"))
                            or ($labels["goog-batch-worker"] == "true")
                    )
                ] | length)
            ] | @tsv
        ' <<< "$current_json"
    )"
    IFS=$'\t' read -r total_count new_count batch_count <<< "$counts"

    if [[ -t 1 ]]
    then
        printf '\033[2J\033[H'
    fi
    printf '%s  total=%s  new_since_start=%s/%s  batch_workers=%s\n' \
        "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" \
        "$total_count" "$new_count" "$max_new_vms" "$batch_count"
    printf 'STATE\tNAME\tZONE\tMACHINE_TYPE\tBATCH_JOB\tCREATED\n'
    jq -r --argjson baseline "$baseline_names" '
        .[]
        | .name as $name
        | [
            (if ($baseline | index($name)) == null then "NEW" else "BASELINE" end),
            .name,
            (.zone | split("/")[-1]),
            (.machineType | split("/")[-1]),
            ((.labels // {})["batch-job-id"] // "-"),
            (.creationTimestamp // "-")
        ]
        | @tsv
    ' <<< "$current_json"

    if ((new_count > max_new_vms))
    then
        printf '\aALERT: %s new running VMs exceeds the expected maximum of %s.\n' \
            "$new_count" "$max_new_vms" >&2
    else
        printf 'OK: new running VM count is within the expected maximum.\n'
    fi

    sleep "$interval_seconds"
done
