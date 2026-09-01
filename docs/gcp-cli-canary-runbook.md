# GCP CLI canary controller

This runbook prepares a clean controller for a bounded RNA-seq cloud canary.
It does not use the OmicsPipelines GUI, and workflow submission remains an
explicit operator-gated step.
The controller authenticates as the existing `cromwell-prod` service account;
the operator's local account does not need direct Google Batch permissions.

## 1. Start the VM guard

Run this from the modernization repository in a dedicated local terminal before
creating the controller:

```bash
bash scripts/gcp/watch_running_vms.sh motrpac-portal 4 15
```

The baseline is whatever is already running when the watcher starts. The limit
of four new VMs is the benchmark cap: one controller plus at most three
concurrent Batch workers. Launch with the reviewed benchmark configuration's
three-job concurrency limit. Any unrelated VM started later is conservatively
counted as new. Batch workers are also identified by Google's `batch-node` and
`batch-job-id` labels and its `goog-batch-worker` marker.

## 2. Create the inert controller

Run this command yourself in a second local terminal. It creates one standard,
on-demand controller with no startup script; creating it cannot submit a Batch
job.

Do not repurpose `motrpac-data-file-workload` for this role. It is shared
data-deposition infrastructure with a different service account and lifecycle.

```bash
gcloud compute instances create rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --machine-type=e2-standard-4 \
  --provisioning-model=STANDARD \
  --network=default \
  --subnet=default \
  --service-account=cromwell-prod@motrpac-portal.iam.gserviceaccount.com \
  --scopes=cloud-platform \
  --image=ubuntu-2204-jammy-v20260826 \
  --image-project=ubuntu-os-cloud \
  --boot-disk-size=50GB \
  --boot-disk-type=pd-balanced \
  --boot-disk-auto-delete \
  --no-restart-on-failure \
  --no-deletion-protection \
  --labels=purpose=rnaseq-modernization,stage=benchmark
```

The watcher should now show exactly one new non-Batch VM.

## 3. Copy and verify the reviewed revision

Transfer a complete Git bundle rather than recursively copying the checkout.
This works for ordinary clones and linked worktrees, and transfers only the
committed, reviewable history.

```bash
cd ~/repos/rnaseq-pipeline-audit/worktrees/modernization-production
test -z "$(git status --porcelain)"
git bundle create /tmp/motrpac-rnaseq-modernization-production.bundle modernization-production
git bundle verify /tmp/motrpac-rnaseq-modernization-production.bundle
```

```bash
gcloud compute scp \
  /tmp/motrpac-rnaseq-modernization-production.bundle \
  rnaseq-modernization-controller:~/ \
  --project=motrpac-portal \
  --zone=us-west1-a
```

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='git clone --branch modernization-production ~/motrpac-rnaseq-modernization-production.bundle ~/motrpac-rna-seq-pipeline-modernization && git -C ~/motrpac-rna-seq-pipeline-modernization status --short --branch && git -C ~/motrpac-rna-seq-pipeline-modernization rev-parse HEAD'
```

## 4. Install only the controller prerequisites

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='sudo apt-get update && sudo apt-get install -y ca-certificates curl jq openjdk-17-jre-headless python3-venv'
```

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='mkdir -p ~/tools && cd ~/tools && curl -fL --retry 3 -o cromwell-92.jar https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar && echo "e0e3a050d4124e81369a79059e5774142b2f06bd89df4a0b035f559db85cedf5  cromwell-92.jar" | sha256sum -c - && java -version'
```

Install the repository's existing input-generator requirements in an isolated
controller environment:

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='python3 -m venv ~/tools/rnaseq-controller && ~/tools/rnaseq-controller/bin/pip install -r ~/motrpac-rna-seq-pipeline-modernization/scripts/requirements.txt && ~/tools/rnaseq-controller/bin/python -c "import gcsfs; print(gcsfs.__version__)"'
```

## 5. Verify the attached identity and Batch visibility

The first command must print
`cromwell-prod@motrpac-portal.iam.gserviceaccount.com`:

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='curl -fsS -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/email && echo'
```

The second command performs a read-only Batch API query using that attached
identity. It does not require `batch.*` permissions on the local user account:

```bash
gcloud compute ssh rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a \
  --command='token=$(curl -fsS -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token | jq -r .access_token); curl -fsS -H "Authorization: Bearer ${token}" "https://batch.googleapis.com/v1/projects/motrpac-portal/locations/us-west1/jobs?pageSize=1" | jq "{visible_jobs: ((.jobs // []) | length)}"'
```

Stop the controller without deleting its disk while the reviewed workflow input
is being prepared:

```bash
gcloud compute instances stop rnaseq-modernization-controller \
  --project=motrpac-portal \
  --zone=us-west1-a
```

## 6. Select published assets and an explicit runtime profile

The published asset catalog is
`config/backends/gcp/benchmark-assets-v1.json`. It records the monitoring
script and the deterministic 100k and 5M paired-input sets with their sizes,
SHA-256 values, and GCS object generations. The v47 release profile now points
to checksum-addressed references and digest-addressed images.

Runtime sizing is deliberately explicit; the generator does not infer a
profile from a filename or compressed byte count:

| Profile | Reviewed envelope | STAR CPU/RAM/disk GB/type | RSEM CPU/RAM/disk GB |
| --- | --- | --- | --- |
| `runtime-human-v47-small-v1.json` | One sample, at most 5M pairs | 10 / 64 / 90 / HDD | 10 / 16 / 30 |
| `runtime-human-v47-full-lean-v1.json` | Accepted v47 canary and full-sample benchmark envelope | 12 / 72 / 120 / SSD | 10 / 40 / 60 |
| `runtime-human-v47-high-candidate-v1.json` | Manual candidate derived from one 54.4M-pair, 7.92-GB adipose canary | 12 / 72 / 150 / SSD | 10 / 40 / 60 |

These are buffered starting points from local and cloud profiling evidence,
not final cohort-wide recommendations. The lean profile and 120 GB SSD passed
both the 5M canary and the tested 39.4-million-pair full sample. The generator
continues to select HDD when `--star-disk-type` is omitted, preserving historical
behavior; benchmark inputs must pass `--star-disk-type SSD` explicitly. STAR
retains substantial memory for its fixed index footprint, and Picard
MarkDuplicates retains 36 GB in the lean profile. Heterogeneous production
samples are still needed before fitting a size-dependent provisioning rule.
The high-input candidate retains every lean CPU and memory value and raises
only STAR scratch. Its 150-GB size is inferred from a 110.63-GiB sampled peak
on the adipose canary and is projected to leave about 35 GiB free. The exact
candidate still requires one production-interface canary before broader use.

The project's `default` VPC uses custom subnet creation. The Batch backend
therefore pins the `default` network and resolves its `default` subnet in each
task's worker region; omitting that subnet prevents worker creation. The Batch
API location does not select the worker region. With no `zones` override, the
Cromwell default remains `us-central1-b`, matching historical production jobs.
Workers retain external IP addresses for the public Quay images in this first
benchmark.

## 7. Generate reviewed canary inputs without submitting them

Generate the 100k input from the repository root on the controller:

```bash
mkdir -p /tmp/rnaseq-gcp-100k
~/tools/rnaseq-controller/bin/python scripts/make_json_rnaseq.py \
  --gcp_path gs://omicspipelines-get/benchmarks/rnaseq-modernization/v1/inputs/13146031002/source-order-v1/100k \
  --output_path /tmp/rnaseq-gcp-100k \
  --output_report_name muscle_13146031002_100k_v47 \
  --organism human \
  --version gencode_v47 \
  --num_chunks 1 \
  --index \
  --combine-contamination-qc \
  --contamination-qc-pairs 100000 \
  --runtime-profile config/backends/gcp/runtime-human-v47-full-lean-v1.json \
  --star-disk-type SSD \
  --project motrpac-portal
```

Use the same command for 5M after changing `100k` to `5m` in both paths,
changing the report name, and setting `--contamination-qc-pairs 1000000`.
The full muscle input uses its recorded raw-data prefix and the same lean
profile and SSD selection.

For the on-demand procedure in this runbook, before submission compare every
controlled-input generation with the asset
catalog, inspect the generated JSON, confirm that
`rnaseq_pipeline.num_preemptible_attempts` is absent or zero, restart the VM
watcher, review the exact Cromwell command, and record the clean launch commit:

```bash
git rev-parse HEAD > ~/rnaseq-canary-launch-revision.txt
```

Do not submit a workflow merely as a consequence of generating these inputs.
The one-off Spot-to-on-demand
launcher is retired rather than part of the production path. Any new controlled
benchmark package must be checksum-reviewed, obtain its workflow ID from the
structured `cromwell-workflow-id` label on the first Batch job, validate the
top-level attempt preemptibility plus Batch provisioning model, and accept its
reviewed run name before launch. The normal Caper/server path should instead use
a unique Caper string label and Cromwell's JSON workflow-query API. Neither path
should parse human log text.

## 8. Capture evidence after completion

Export the workflow's complete Cromwell metadata JSON through the normal
operator interface, then run the capture on the controller. Record the clean
checkout revision at launch as shown above; capture refuses
to attribute the run to a different checkout. Do not run it under the local
user account: that account cannot describe Batch jobs, while the controller's
attached `cromwell-prod` identity can.

```bash
cd ~/motrpac-rna-seq-pipeline-modernization
test -z "$(git status --porcelain)"
bash scripts/gcp/capture_workflow_evidence.sh \
  /path/to/complete-metadata.json \
  /path/to/new-evidence-directory \
  "$(cat ~/rnaseq-canary-launch-revision.txt)"
jq -e '
  .workflow_status == "Succeeded" and
  .complete == true and
  .top_level_output_object_count == 12
' /path/to/new-evidence-directory/capture-status.json
python3 scripts/gcp/summarize_workflow_cost.py \
  /path/to/new-evidence-directory \
  --output /path/to/new-cost-summary.json
```

The checkout must remain clean from submission through capture; the evidence
bundle records both its exact revision and cleanliness state.

The capture fails if any attempt, Batch job, monitoring/stdout/stderr stream,
or immutable input/output object metadata is missing. It copies only the
explicitly allowlisted top-level matrices, QC report, contamination manifest,
and UMI manifests/metrics, with a 256-MiB per-object ceiling; it never walks an
execution root. `evidence-manifest.sha256` covers every captured file. The cost
summary validates that manifest and rolls attempts into the fixed eight
scientific pipeline phases as well as Batch lifecycle phases. `complete` means
that the listed transport evidence was captured; the separate command above
requires the successful full-canary 12-object output contract. Apply the
existing matrix and QC scientific gates before acceptance. Transfer the
completed directory and cost summary before stopping the controller.

## Handoff gate

The authorized next scope is one clean, human-v47, on-demand, full-depth canary
through the operator's normal interface. Before submission, confirm that the
interface is using this checked GCP Batch backend and monitoring configuration;
stop if it still targets PAPI or another backend. Name the sample and profile in
the launch record: use the lean profile only within its tested 39.4-million-pair,
5.79-GB envelope, or explicitly review the high-input candidate. Disable call
cache/reference disks and allow at most three Batch workers. Do not release the
remaining cohort until its outputs, evidence bundle, cost summary, telemetry,
and worker cleanup have been reviewed.

General release remains blocked on:

- verification of the real Caper/GUI backend, output ingestion, and monitoring;
- accepted modern image/reference profiles for v39 and rat; and
- a production Spot-only/market-policy implementation and bounded retry test.
