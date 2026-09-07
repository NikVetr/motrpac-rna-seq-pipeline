#!/usr/bin/env bash
# Run locally after copying evidence off a dedicated, single-workflow controller.
set -euo pipefail

if [[ "$#" -ne 4 ]]; then
    echo "usage: $0 LOCAL_EVIDENCE_DIRECTORY PROJECT DEDICATED_CONTROLLER ZONE" >&2
    exit 2
fi

evidence="$1"
project="$2"
controller="$3"
zone="$4"

jq -e '.complete == true and .workflow_status == "Succeeded" and .missing_artifact_count == 0' \
    "$evidence/capture-status.json" >/dev/null
(cd "$evidence" && sha256sum --check --quiet evidence-manifest.sha256)

gcloud compute instances stop "$controller" --project="$project" --zone="$zone" --quiet
status="$(gcloud compute instances describe "$controller" --project="$project" \
    --zone="$zone" --format='value(status)')"
[[ "$status" == TERMINATED ]] || { echo "controller is not TERMINATED: $status" >&2; exit 1; }
