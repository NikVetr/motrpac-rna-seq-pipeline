# Task profiling for cohort calibration

The monitoring asset records cgroup anonymous/file/inactive-file memory alongside
CPU, total memory and disk use. It is the immutable script used in the v47/v50
cloud pilot. `summarize_workflow_cost.py` reports peak working memory as the
maximum of `current - inactive_file` at each sample, preserves cache-inclusive
peaks, and accepts older logs. Working memory is not minimum required RAM.

Use the `monitoring_script` URI from
`config/backends/gcp/workflow-options-benchmark.example.json` in actual submission
options. The worker account needs read access to that regional bucket. Keep
worker/execution locality aligned. The example's disabled caching and bounded
backend settings are not production defaults: enable persistent caching and
`ContinueWhilePossible` for calibration repairs.

Export full Cromwell metadata including call inputs, outputs and runtime
attributes, then capture a snapshot into a new directory:

```bash
python3 scripts/gcp/capture_task_profiles.py metadata.json profiles/run-001
```

The standard-library collector uses gcloud Storage object-read and Batch
job-read permissions. It launches no workers. It describes input/output objects,
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
produce an explicit incomplete-capture error. Live/aborted streams may never
have been exported; unavailable measurements must not be treated as final peaks
or omitted from failure accounting. Capture before deleting jobs or execution
files. Use a separate directory and metadata export for every snapshot.

The existing complete-workflow evidence/cost tools retain their documented
scope. This collector supports larger and partial workflows without downloading
their final matrices. Join its attempt identities with actual billed durations
and market rates when preparing cohort cost summaries.
For us-west2 N1 jobs, the matching frozen rates are in
`config/backends/gcp/gcp-rates-n1-us-west2-20260911.json`; pass this through
`summarize_workflow_cost.py --rates` rather than using its Americas default.

Start with two full-depth v50 libraries within the intended approximately
100-library set, check output/profiling collection, then expand while retaining
those results. Include prior UMI/RSEM extremes and missing-UMI samples where
applicable. Another paired v47 run is unnecessary. The smaller RAM candidates
still need cloud validation; the monitor has already run in the cloud.

Compare actual inputs against working/anonymous memory, scratch, duration and
allocated CPU/RAM. Include failures and retries. Compare RSEM retained templates
and alignment records alongside BAM bytes, UMI incoming work, and each task's
localized bytes. Retain reference and expression-mode identity. CPU-thread
changes require separate runtime/cost measurements; this pilot did not test
speedup from additional threads.
