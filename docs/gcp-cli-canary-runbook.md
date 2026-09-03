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

The supplied backend accepts at most 10 samples in a scatter and 200 task calls
in its single active workflow. Those are graph-size safeguards, not concurrency
settings. The separate three-job limit remains the hard worker ceiling; the
larger graph limits let a 5--10-sample pilot finish instead of rejecting it at
the original one-sample canary caps.

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

## 6. Select published assets and a runtime floor

The published asset catalog is
`config/backends/gcp/benchmark-assets-v1.json`. It records the monitoring
script and the deterministic 100k and 5M paired-input sets with their sizes,
SHA-256 values, and GCS object generations. The v47 release profile now points
to checksum-addressed references and digest-addressed images.

Runtime profiles set fixed CPU and memory requests and a minimum STAR scratch
allocation:

| Profile | Intended use | STAR CPU/RAM/minimum disk GB/type | RSEM CPU/RAM/disk GB |
| --- | --- | --- | --- |
| `runtime-human-v47-small-v1.json` | One sample, at most 5M pairs | 10 / 64 / 90 / HDD | 10 / 16 / 30 |
| `runtime-human-v47-full-lean-v1.json` | Full-depth human v47 samples | 12 / 72 / 120 / SSD | 10 / 40 / 60 |
| `runtime-human-v47-high-candidate-v1.json` | Conservative fixed 150-GB minimum | 12 / 72 / 150 / SSD | 10 / 40 / 60 |

After Cutadapt, the workflow reads the exact number of surviving pairs for each
sample and raises that sample's STAR disk independently:

| Post-trim read pairs | STAR scratch GB |
| ---: | ---: |
| at most 5 million | 90 |
| at most 40 million | 120 |
| at most 65 million | 150 |
| at most 90 million | 180 |
| at most 110 million | 200 |
| at most 155 million | 250 |
| at most 200 million | 300 |
| more than 200 million | 400 |

The effective allocation is the larger of this tier and the profile's
`star_disk` value. Thus one input JSON may contain heterogeneous samples; the
operator does not need to generate a separate JSON for each disk tier. The
120-GB tier is retained for samples of at most 40 million post-trim pairs. Its
direct storage saving relative to 150 GB is only about $0.005 for a 41-minute
STAR call at the reviewed `us-west1` SSD rate, so it is useful but not a major
cost lever.

The lean profile and 120-GB SSD passed the tested 39.4-million-pair full sample
and the operator-interface canary. A separate 150-GB profile passed the
48.1-million-pair CLI canary. The tiers add conservative headroom above the
observed relationship and return to the historical 400-GB allocation above
200 million pairs. CPU and memory remain fixed because the v47 index establishes
a large input-independent floor and the current evidence does not support a
more complicated policy. The generator continues to select HDD when
`--star-disk-type` is omitted, preserving historical behavior; production v47
inputs should pass `--star-disk-type SSD`.

The project's `default` VPC uses custom subnet creation. The Batch backend
therefore pins the `default` network and resolves its `default` subnet in each
task's worker region; omitting that subnet prevents worker creation. For the
next pilot, set both the GCP Batch job location and worker-zone preference to
`us-west1`; current GCP Batch location rules require the allowed worker zones
to be in the job's region. Workers retain external IP addresses for the public
Quay images in this benchmark phase.

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
  --run-multiqc \
  --combine-contamination-qc \
  --contamination-qc-pairs 100000 \
  --runtime-profile config/backends/gcp/runtime-human-v47-full-lean-v1.json \
  --star-disk-type SSD \
  --project motrpac-portal
```

Use the same command for 5M after changing `100k` to `5m` in both paths,
changing the report name, and setting `--contamination-qc-pairs 1000000`.
For a bounded multi-sample pilot, create a local manifest containing the exact
FASTQ prefixes, one per line, and add:

```bash
  --sample-list /path/to/pilot-samples.txt \
  --num_chunks 1
```

The manifest selects only those samples from the single GCS prefix passed to
`--gcp_path`. One JSON contains the full pilot, and each scattered STAR call
receives its own post-Cutadapt disk tier. After acceptance, reuse the same
manifest with `--exclude-sample-list` to generate a nonoverlapping remainder.
Write the remainder JSONs to a fresh output directory; the generator refuses to
mix a new generation with existing `set*_rnaseq.json` files. Do not combine the
two selection flags.

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
  .missing_artifact_count == 0
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
UMI manifests/metrics, and optional MultiQC archives, with a 256-MiB per-object
ceiling; it never walks an execution root. `evidence-manifest.sha256` covers
every captured file. The cost summary validates that manifest and rolls
attempts into the fixed eight scientific pipeline phases as well as Batch
lifecycle phases. `complete` means that the listed transport evidence was
captured. A multi-sample workflow has a variable output-object count because
per-sample QC and provenance arrays grow with the scatter. Apply the existing
matrix and QC scientific gates before acceptance.
Transfer the completed directory and cost summary before stopping the
controller.

## Handoff gate

The authorized next scope is one clean, human-v47, on-demand pilot of 5--10
pre-/3.5-hour muscle samples through the operator's normal interface. Select
paired timepoints spanning the available compressed-input sizes from one
confirmed GCS prefix, record them in an exact sample manifest, and use the lean
profile with SSD scratch. Before submission, confirm that the interface uses
the checked GCP Batch backend and monitoring configuration in `us-west1`; stop
if it targets PAPI or another backend. Enable MultiQC, disable call
cache/reference disks, retain the reviewed worker-concurrency limit, and pass
`--combine-contamination-qc --contamination-qc-pairs 1000000`. Sampled
contamination QC remains a pilot setting until its estimates have been checked
against the corresponding historical full-depth results. Do not release the
remaining cohort until canonical UMI outputs, the default absence of all-read
outputs, dynamic disk tiers, evidence capture, cost summary, telemetry, and
worker cleanup have been reviewed. Reuse the pilot manifest as an exclusion
list for the subsequent cohort input.

General release remains blocked on:

- multi-sample verification of canonical-output ingestion and dynamic sizing;
- accepted modern image/reference profiles for v39 and rat; and
- a production Spot-only/market-policy implementation and bounded retry test.
