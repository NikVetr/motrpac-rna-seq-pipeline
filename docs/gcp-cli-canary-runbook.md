# GCP CLI canary

Use a dedicated controller and a small, explicit sample set for isolated
resource measurements. For production cohorts and app integration, use
[cohort provisioning](cohort-provisioning.md).

## Controller and backend

Use an existing controller with Java 17, Python 3, gcloud, jq and the input
generator dependencies in `scripts/requirements.txt`. Install the Cromwell
release in `config/backends/gcp/cromwell-release.json` and verify its SHA-256.
Check out the intended pipeline commit and record `git rev-parse HEAD` at
submission. Capture evidence from that same revision.

Controller credentials need Batch job submission/read access in the target
project, access to the worker service account, and the required GCS permissions.
Workers need read access to inputs, references, images and the monitoring
script, and write access to execution storage. Configure the network/subnet for
the target project and region.

The supplied `config/backends/gcp/google_batch.conf` runs one workflow with at
most 10 samples, 200 task calls, three concurrent workers and a four-hour task
timeout. It disables caching for cold measurements. Adapt these explicit limits
to the intended experiment; do not use this backend for a full cohort.

The backend reads these environment variables:

| Variable | Value |
|---|---|
| `MOTRPAC_GCP_PROJECT` | Worker project |
| `MOTRPAC_GCP_COMPUTE_SERVICE_ACCOUNT` | Worker service-account email |
| `MOTRPAC_GCP_BATCH_ROOT` | Regional GCS execution prefix |
| `MOTRPAC_GCP_BATCH_LOCATION` | Matching worker region |

The read-only VM watcher reports new VMs relative to its starting baseline:

```bash
bash scripts/gcp/watch_running_vms.sh PROJECT 4 15
```

Four allows one new controller plus three workers. Unrelated new VMs also count;
the watcher reports limits but does not stop resources.

## Inputs and launch

`config/backends/gcp/benchmark-assets-v1.json` records controlled 100k/5M
input sets and monitoring-script checksums/generations. Verify asset availability
and locality before reuse. For full-depth samples, select the appropriate
human or rat profile and an exact sample list; do not use the small v47 profile
for full-depth work.

```bash
python3 scripts/make_json_rnaseq.py \
  -g gs://INPUT_BUCKET/BATCH/fastq_raw -o input_json -r canary \
  -a human -v gencode_v50 -n 1 -p PROJECT \
  --sample-list selected_samples.txt \
  --runtime-profile config/backends/gcp/runtime-human-v50-full-candidate-v1.json \
  --star-disk-type SSD --umi-dup-disk-type SSD \
  --combine-contamination-qc --contamination-qc-pairs 1000000 -i
```

Copy `config/backends/gcp/workflow-options-benchmark.example.json` to
`options.json`; set worker zones and a readable monitoring-script URI for the
chosen region. Cold options disable cache reads/writes. The workflow defaults
to one Spot attempt followed by on-demand; set
`rnaseq_pipeline.num_preemptible_attempts` to `0` in the generated input JSON
for an on-demand comparison.

With the backend environment variables set, run:

```bash
python3 scripts/gcp/check_locality.py \
  --inputs input_json/set1_rnaseq.json --options options.json -- \
  java -Dconfig.file=config/backends/gcp/google_batch.conf -jar cromwell-92.jar \
  run wdl/rnaseq_pipeline_scatter.wdl \
  -i input_json/set1_rnaseq.json -o options.json
```

Omit the command after `--` for a read-only locality check. The guard rejects
execution/worker mismatches and reports remote sources; it does not move data.
For Caper, use its server configuration and a unique submission label instead
of this standalone launch.

## Evidence and cleanup

Export complete Cromwell metadata, retaining inputs, outputs and all attempts.
Use [task profiling](task-profiling.md) for large or incomplete cohorts.
For a complete canary, use an identity with Batch job-read and GCS permissions:

```bash
bash scripts/gcp/capture_workflow_evidence.sh \
  metadata.json evidence SUBMITTED_COMMIT
jq -e '.workflow_status == "Succeeded" and .complete == true and .missing_artifact_count == 0' \
  evidence/capture-status.json
```

The capture checks the submission revision, records checkout cleanliness,
describes input/output objects, and captures Batch jobs and task logs. It
downloads only allowlisted small outputs with a 256-MiB per-object limit.
A SHA-256 manifest covers captured files. Complete capture establishes evidence
availability; separately verify matrix values, sample order, QC and expression
metadata before accepting outputs.

`scripts/gcp/summarize_workflow_cost.py` consumes this evidence with an explicit
matching `--rates` manifest. Its current rate formats support N1/N2 and listed
predefined N1 shapes, not E2. For E2, retain task profiles and actual allocations
for billing/rate analysis; do not apply N1 rates.

Copy accepted outputs and evidence before stopping a dedicated controller:

```bash
bash scripts/gcp/stop_controller_after_capture.sh LOCAL_EVIDENCE PROJECT CONTROLLER ZONE
```

Do not use this helper on a shared Caper server. Retained disks and objects
continue to incur storage charges; review intermediate retention separately.
