# Cohort provisioning and operations

For the per-attempt input, count and memory evidence needed to recalibrate
these rules, see [task profiling](task-profiling.md).
For version-specific human v50 settings and cache-enabled submission options,
see the [initial v50 cohort calibration](v50-cohort-calibration.md).

Resource checks use all 4,455 successful per-sample monitoring logs from the
296-library workflow `68aaf86f-b391-4423-a269-21371a91f416` and single-library
workflow `6891b68f-1751-4c8f-8acb-b5049377fa9d`. Logs were refreshed using final
Cromwell attempt identities, with GCS generations recorded. Partial snapshots
must not be used for resource ceilings or final retry-cost accounting.

## Resource rules

All configured resource requests remain floors; CPU counts are unchanged.
The table below describes the established v47 full-depth profile.

| Task | Rule | Evidence and scope |
|---|---|---|
| STAR | Existing post-trim read-count scratch tiers | 297 successful libraries; peak memory 62.06 GiB. Keep the 72-GB full-depth RAM setting. |
| UMI | Existing `ceil(2 * input_BAM_GiB + 15)` scratch with 80-GB full-depth floor | Peak memory 18.77 GiB. Keep 20 GB RAM and the SSD floor; the cohort does not establish SSD performance below that floor. |
| Molecule RSEM | RAM `4 * ceil((16 + 2 * input_BAM_GiB) / 4)`; scratch `ceil(10 + 4 * input_BAM_GiB)` | Peak memory 38.26 GiB. With 40-GB RAM and 60-GB disk floors, the rules cover all 297 samples with at least 5.82 GiB memory and 15.62 GiB disk headroom. Only three samples request more than 40 GB RAM (44–48 GB). |
| Both merges | RAM `4 * ceil(libraries / 75)`; existing scratch `ceil(3 * input_GiB + 10)` | The successful 297-library merge used 16 GB RAM and 68 GB disk; observed peaks were 9.94 and 38.39 GiB. |

Merge input directories contain symbolic links to localized files. The merge
helpers and image are unchanged; no second input copy is needed. Keep the
validated disk allowance until a cloud run establishes the new disk peak.

RSEM growth and merge RAM rules are conservative capacity rules calibrated on
human v47, not universal bounds for other annotations or much larger cohorts.
The small profile remains a small-run floor; it does not reduce the full-depth
profile. Do not reduce configured floors from this retrospective analysis.
The WDL expressions and link-based merge commands pass local checks and a
three-call GCP canary (`ef97bfae-cc72-40cf-9d5d-869ea4134c39`). One RSEM rerun
used 4 vCPUs, 20 GB RAM, and 15 GB scratch; both two-library merges used
1 vCPU, 4 GB RAM, and 11 GB scratch. All outputs matched their original values
exactly, with one worker at a time and no retries. The canary used 1.68 worker
vCPU-hours and $0.0957 frozen-rate modeled worker cost, excluding controller,
transfer, and retained storage. Large-sample performance and smaller SSDs
remain unbenchmarked; no CPU speedup is claimed from average utilization.

## Locality belongs at launch

Use a regional execution bucket in the same region as Batch and all permitted
worker zones. This prevents repeated intermediate transfers between regions.
Upstream demultiplexed FASTQs can remain in their authoritative bucket.
Stage references or FASTQs only when repeated consumption justifies a second
copy; preserve object identity and update input URIs explicitly.

The operator execution bucket `omicspipelines-test-data` is in `us-west2`.
The recovery bucket `omicspipelines-get` is `US` multi-region and is not an
explicitly colocated regional execution bucket. A different project's worker
service account also needs permission to use the chosen bucket.

Gate a local Cromwell launch using the same inputs/options passed to Cromwell:

```bash
env MOTRPAC_GCP_BATCH_ROOT=gs://REGIONAL_BUCKET/executions \
    MOTRPAC_GCP_BATCH_LOCATION=us-west2 \
    python3 scripts/gcp/check_locality.py \
    --inputs input.json --options options.json -- \
    java -Dconfig.file=config/backends/gcp/google_batch.conf -jar cromwell.jar \
    run wdl/rnaseq_pipeline_scatter.wdl -i input.json -o options.json
```

Set `default_runtime_attributes.zones` in those options to
`us-west2-a us-west2-b us-west2-c`. Existing project/service-account variables
and cohort concurrency settings must also be supplied. The checked-in backend
is a bounded benchmark configuration, not a 297-library launch profile.

Set `"workflow_failure_mode": "ContinueWhilePossible"` in submission options,
as in the benchmark example, so independent sample tasks continue after a
terminal task failure. This is a submission option and needs no server restart.
The workflow still reports failure, and a merge requiring the failed sample
waits for repair. Preserve execution outputs and enable persistent call caching
for production repair runs; this option does not retry failures with more RAM.

The guard rejects mismatched execution storage, zones, and workflow root
overrides; it reports remote input/reference buckets without moving them.
Running it without a trailing command performs a read-only check. With Caper
server submissions, configure the server's backend root and region first:
environment variables on the submitting client do not reconfigure that server.
Direct engine launches can bypass this guard, so use the guarded entry point
in launch automation. The WDL itself remains portable and cloud-independent.

## Metadata and retention

Prepare downstream metadata from explicit study records and modern merged QC:

```bash
python3 scripts/prepare_sample_metadata.py \
    --matrix rsem_genes_count.txt --study study_metadata.csv --qc cohort_qc.csv \
    --output sample_metadata.csv \
    --required-study-columns pid visitcode Timepoint randomGroupCode Sex \
    calculatedAge BMI codedsiteid Batch RIN
```

Study records may contain other samples; QC must match the matrix sample set
exactly. Duplicate IDs, overlapping old QC columns, missing required values,
and existing outputs fail. Values remain strings and are not rescaled.

Keep final matrices, sample/study provenance, manifests, QC, and complete
checksummed evidence. Review bulky execution intermediates only after output
acceptance and any planned comparisons; deleting them invalidates cache and
recovery references. Apply any lifecycle rule to a dedicated disposable prefix,
never broadly to raw inputs or reference buckets. No deletion policy is applied
automatically by these changes.

For a dedicated single-workflow controller, copy evidence locally first, then:

```bash
bash scripts/gcp/stop_controller_after_capture.sh LOCAL_EVIDENCE PROJECT CONTROLLER ZONE
```

This verifies successful, complete, checksummed evidence before stopping the
explicit controller and confirming `TERMINATED`. Do not use it on a shared
Caper server. Retained disks and objects continue to incur storage charges.

## Checks

Run `python3 -m unittest discover -s tests`, then, in the existing MiniWDL
environment, `python tests/check_wdl_resources.py`. The latter evaluates the
actual WDL resource expressions and executes both rendered merge commands
on fixtures, comparing outputs byte for byte. Validate the full graph with
MiniWDL and WOMtool before cloud submission.
