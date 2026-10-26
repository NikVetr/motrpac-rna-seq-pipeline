# Task profiling for cohort calibration

The monitoring asset records cgroup anonymous/file/inactive-file memory alongside
CPU, total memory and disk use. `summarize_workflow_cost.py` reports peak working
memory as the maximum of `current - inactive_file` at each sample, preserves cache-inclusive
peaks, and accepts older logs. Working memory is not minimum required RAM.

Use the `monitoring_script` URI from
`config/backends/gcp/workflow-options-cohort.example.json` in actual submission
options. The worker account needs read access to that regional bucket. Keep
worker/execution locality aligned. The cohort options enable cache reads/writes
and `ContinueWhilePossible`; persistent server-side caching is also required.
The separate cold-benchmark options and bounded backend are not cohort settings.

Export full Cromwell metadata including call inputs, outputs and runtime
attributes, then capture a snapshot into a new directory:

```bash
python3 scripts/gcp/capture_task_profiles.py metadata.json profiles/run-001
```

The standard-library collector uses the active gcloud identity and requires Storage
object-read and Batch job-read permissions. The default REST reader obtains one
access token (refreshing it as needed), reuses HTTPS connections across worker
pools and shares identical reads within the cohort capture. Successful responses
are cached on temporary local disk and discarded when the command exits.
`--workers 8` bounds concurrent reads; `--transport gcloud` explicitly selects
the subprocess reader. There is no automatic permission-error fallback.
The REST reader supports direct HTTPS and HTTP CONNECT proxies; other proxy
schemes require the gcloud transport.
It launches no cloud workers. It describes input/output objects,
captures actual Batch allocations and attempt states, and downloads only logs
and small existing metrics. BAMs, FASTQs, reference archives and expression
matrices are never downloaded. Input parameter names keep reference bytes
separate from BAM/FASTQ bytes. Scalar outputs retain post-trim counts; UMI metrics
retain incoming/retained templates and alignment counts; RSEM `.cnt` statistics
cover molecule and all-read runs. No extra scientific task or BAM scan is needed.

`profiles.json` records terminal attempts, failures, caching and pending-attempt
identities. `objects.json` retains size/generation; downloads use that generation.
A SHA-256 manifest covers captured evidence. Failed parent workflows are
supported without requiring a final merge. Missing objects are recorded and
produce an explicit incomplete-capture error, except absent logs (monitoring
log, stdout, stderr) of cache hits, which never ran, and of attempts that did not
finish (preempted, failed or aborted): these print a warning and keep
`complete_capture: false`. Successful attempts that ran must have their logs.
Workflow files are sized with one listing per workflow root. Each retry retains its own logs; successful-attempt
measurements never replace missing earlier measurements. Batch read failures
also preserve the remaining evidence before reporting an error. Live/aborted streams may never
have been exported; unavailable measurements must not be treated as final peaks
or omitted from failure accounting. Capture before deleting jobs or execution
files. Use a separate directory and metadata export for every snapshot.

For final cohort preservation, capture every run referenced by the merge inputs
and recursively follow its cache origins:

```bash
python3 scripts/gcp/capture_cohort_profiles.py profiles/current-tissue \
  --merge-inputs merge-inputs.json --workers 8 \
  --reuse-captures profiles/previous-tissue
```

Omit `--reuse-captures` for the first capture. It can be repeated for several
previous local cohort directories, including unpacked preserved bundles.
Each reused workflow must match freshly read metadata exactly and pass all
SHA-256 checks, terminal-state and attempt-coverage checks. A full capture also
requires the stored Batch description for every attempt with a job ID.
Changed metadata or missing Batch descriptions cause a fresh capture; corrupt
evidence fails. Reused files are copied into the new self-contained bundle, and
`runs.json` records their source directories. This avoids repeatedly downloading
the same shared calibration workflow across tissues. It does not assert that
old execution objects still exist; cleanup must check their current generations
and dependencies separately.

The cohort command records per-workflow failures in `runs.json`, including missing
files in an unpacked reuse source, and continues checking the remaining workflows.
It refuses running workflows or pending attempts and creates no bundle when a
capture fails. Expected absent cache-hit/interrupted-attempt logs
remain explicit and do not prevent bundling. Gzip level 1 reduces compression
work while retaining the same evidence. `--destination gs://.../task-profiles`
uploads only after capture succeeds and refuses to replace an existing object.
Do not use `--skip-batch` to work around missing permissions for a final cost
evidence capture; run with an identity that has `batch.jobs.get` instead.

Start with at most two cohort processes and eight readers each on the operator
controller. Measure controller load and examine actual errors before increasing
concurrency; this is an operational starting point, not a documented GCP limit.
Let existing captures finish before changing their code. Validate one new bundle
and its checksums before adopting the updated reader for the remaining cohorts.
Storage reads use Google's [paginated object listing](https://docs.cloud.google.com/storage/docs/json_api/v1/objects/list)
and [generation-specific downloads](https://docs.cloud.google.com/storage/docs/json_api/v1/objects/get);
allocations use [Batch job GET](https://docs.cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs/get).

Per-sample `qc_diagnostics.json` files are included when available. They retain
RSEM convergence/short-pair diagnostics, strand consistency, featureCounts
assignment categories and FastQC flags alongside resource measurements.
RSEM transcript/background and affected-gene convergence TSVs are captured too,
with a 64-MiB limit per table. Their summary is under `rsem.final_iteration` in
the QC JSON; header-only tables indicate no flagged components/genes.

The collector supports larger and partial workflows without downloading final
matrices. Join attempt identities with actual billed durations and market rates
when preparing cohort cost summaries. The complete-workflow cost summarizer
supports N1/N2 rate manifests and explicitly listed predefined N1 shapes; E2
allocations are captured by this collector but require separate rate analysis.
For us-west2 N1 custom and eligible predefined jobs, the matching frozen rates are in
`config/backends/gcp/gcp-rates-n1-us-west2-20260911.json`; pass this through
`summarize_workflow_cost.py --rates` rather than using its Americas default.
The manifest prices `n1-highmem-2` explicitly. Unlisted predefined shapes fail
instead of being charged at custom rates; actual Batch machine identity is retained.

Compare actual inputs against working/anonymous memory, scratch, duration and
allocated CPU/RAM. Include failures and retries. Compare RSEM retained templates
and alignment records alongside BAM bytes, UMI incoming work, and each task's
localized bytes. Retain reference and expression-mode identity. CPU-thread
changes require separate runtime/cost measurements; utilization alone does not
establish speedup from additional threads.
