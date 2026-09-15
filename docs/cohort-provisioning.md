# Cohort provisioning and operations

Runtime profiles provide CPU requests and RAM/disk floors. The workflow raises
resource requests from the data entering each task. See
[human v50](v50-cohort-calibration.md), [rat v116](rat-ensembl116.md), and
[task profiling](task-profiling.md) for annotation-specific floors and measurements.

## Resource rules

Each effective request is the larger of the configured floor and the rule below.

| Task | Dynamic rule |
|---|---|
| STAR scratch | Post-trim pair-count tiers: 90/120/150/180/200/250/300/400 GB at up to 5/40/65/90/110/155/200/>200 million pairs. Human v50 multiplies the tier by 1.30, rounded up. |
| UMI scratch | `ceil(2 * combined_STAR_BAM_GiB + 15)` |
| RSEM RAM | `4 * ceil((16 + 2 * transcriptome_BAM_GiB) / 4)` GB |
| RSEM scratch | `ceil(10 + 4 * transcriptome_BAM_GiB)` GB |
| Gene/expression/isoform merge RAM | `4 * ceil(libraries / 75)` GB |
| Merge scratch | `ceil(3 * input_GiB + 10)` GB |

Merge inputs are symbolic links to localized files, avoiding a second input
copy. RSEM growth applies to the BAM entering either molecule or all-read
quantification. These capacity rules are starting allocations, not universal
bounds for other annotations or arbitrarily large libraries.

The v47 full-depth profile keeps 72 GB STAR RAM, 20 GB UMI RAM/80 GB scratch,
and 40 GB RSEM RAM/60 GB scratch. A 297-library calibration observed peak working
RAM of 62.06/18.77/38.26 GiB for STAR/UMI/RSEM; the gene merge peaked at 9.94 GiB
RAM and 38.39 GiB scratch. These measurements support the v47 floors and do not
establish v50 transcript-merge requirements.

## Deployment and locality

Use regional execution storage matching the Batch region and permitted worker
zones. Keep authoritative FASTQs in their source bucket. Stage remote references
or FASTQs when repeated reads justify a regional copy; retain checksums and
update input URIs explicitly. Project boundaries affect permissions and billing
ownership; region alignment governs GCS transfer locality.

Check the actual server settings and submission files:

```bash
python3 scripts/gcp/check_locality.py \
  --execution-root gs://REGIONAL_BUCKET/executions --region us-west1 \
  --inputs input.json --options options.json
```

Set options' `default_runtime_attributes.zones` to zones in that region.
The guard rejects mismatched execution storage, zones and root overrides, and
reports remote inputs/references. It does not reconfigure a Caper server.

Use Cromwell 92 for task-specific E2/N1 machine selection. Set server-side Batch
project, location, root, worker identity and network/subnet explicitly. Keep a
persistent Cromwell database with `call-caching.enabled = true`, and enable
cache reads/writes in workflow options. Set
`"workflow_failure_mode": "ContinueWhilePossible"` so independent samples
continue after a terminal failure. Merges requiring a failed sample wait for
repair; this setting does not increase RAM or retry command failures.

For 100-library calibration, permit scatter width of at least 100 and enough
total task calls, including retries (5,000 is sufficient for the default graph).
Set concurrency from project quota and the intended spend/turnaround; a limit
of 300 permits 300 active tasks rather than 300 full pipelines. A 48-hour
task timeout allows long tails. The separate `google_batch.conf` is a bounded,
cache-disabled canary backend and must not replace the cohort server config.

OmicsPipelines bundles its own WDL and input generator. App integration must
synchronize the workflow/imports, matched release references and image digests,
resource profiles, reference-release identifiers, and UMI policy. Ensure app
submissions also carry the monitoring script and workflow options. Verify the
Cromwell/database upgrade, then run two intended cohort samples through the app
before expanding. Keep references and images available independently of
disposable execution data.

## Outputs, metadata and retention

Retain canonical raw RSEM gene/isoform results, matrices, QC, expression metadata,
input/options manifests, source revision, and profiling evidence. Join study
covariates to expression-column order using:

```bash
python3 scripts/prepare_sample_metadata.py \
  --matrix rsem_genes_count.txt --study study_metadata.csv --qc cohort_qc.csv \
  --output sample_metadata.csv \
  --required-study-columns pid visitcode Timepoint randomGroupCode Sex \
  calculatedAge BMI codedsiteid Batch RIN
```

Study records may contain other samples; QC must match the matrix sample set.
Duplicate IDs, overlapping QC columns, missing required values and existing
outputs fail. Values remain strings and are not rescaled.

Accepted results need not be rerun solely because CPU/RAM/disk allocations
change. Reuse them through verified cache hits or exclude completed samples and
gather their retained outputs. Retain failed-attempt evidence before repair.
Delete bulky intermediates only after acceptance and planned comparisons;
deletion invalidates cache/recovery references. Keep raw inputs and reference
assets outside execution-data lifecycle rules.
