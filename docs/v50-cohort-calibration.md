# Initial human v50 cohort calibration

Use `config/backends/gcp/runtime-human-v50-full-candidate-v1.json` for the
initial approximately 100 libraries. It preserves the v47 full-depth profile's
threads and disk floors, changes UMI/RSEM/RNA-QC RAM floors to 36/32/12 GB,
and keeps STAR at 12 threads/72 GB. Shared BAM-size growth remains active.
The workflow multiplies STAR's post-trim read-pair scratch tiers by 1.30 only
when `reference_release` is `gencode_v50`, rounding up and respecting larger
explicit floors. v47 and rat resource rules are unchanged.

The profile enables `rnaseq_pipeline.prefer_predefined_n1`. In us-west2, a
2-vCPU/12-GB RNA-QC request uses `n1-highmem-2` (2 vCPUs/13 GiB), approximately
1.2% cheaper in both Spot and on-demand compute than N1 custom. The heap stays
at 12 GB. Other CPU/RAM requests keep the backend's normal sizing. The locality
guard rejects other regions and explicit CPU-platform/machine overrides while
this policy is enabled; disable it when those settings are needed. Use the
pinned Cromwell 92 backend, which supports conditional `gcp` runtime objects.

These are buffered pilot candidates, not validated resource minima. Four
completed v50 RSEM calls used 7.38/27.82/31.85/35.71 GiB working memory and
would request 32/48/52/56 GB. Incomplete UMI tails reached 30.94 GiB working
memory; completed RNA-QC reached approximately 10 GiB. Larger v50 tails and
the cohort transcript merge still need measurements at these allocations.

## Inputs and submission

Generate each input batch using the existing FASTQ location and an explicit
sample list. Retain the scientific settings used for the paired pilot:

```bash
python3 scripts/make_json_rnaseq.py \
  -g gs://BUCKET/BATCH/fastq_raw -o input_json -r v50_calibration \
  -a human -v gencode_v50 -n 1 -p PROJECT \
  --sample-list selected_samples.txt \
  --runtime-profile config/backends/gcp/runtime-human-v50-full-candidate-v1.json \
  --star-disk-type SSD --umi-dup-disk-type SSD \
  --combine-contamination-qc --contamination-qc-pairs 1000000 -i
```

The sample list applies within the specified FASTQ directory; handle multiple
batches with their existing input assembly. Add `--allow-missing-umis` for a
mixed-I1 cohort. Missing-I1 samples remain explicitly marked as not deduplicated.
Default counting is directional UMI molecule expression, with raw gene/isoform
results and matrices retained. A second all-read pass requires
`--retain-all-read-expression`; leave it off for this calibration.

Submit with `config/backends/gcp/workflow-options-cohort.example.json`.
It enables cache reads/writes, `ContinueWhilePossible`, the pinned memory
monitor, us-west2 worker zones, and `maxRetries: 0` for command failures.
The workflow defaults to one Spot attempt, followed by on-demand after a
preemption. In existing JSONs, set `rnaseq_pipeline.num_preemptible_attempts`
explicitly to `1`: an old explicit `0` overrides the new default. Use `2` for
two Spot attempts or `0` for on-demand only. `maxRetries` is independent of
Spot preemption retries; it does not increase RAM after an OOM.

Keep the operator's established regional alignment. Before submission, verify
the actual Caper server has persistent database-backed call caching enabled,
48-hour Batch task timeout, and limits that permit the cohort: scatter width
at least 100, concurrent task VMs initially 20–30, and total jobs per root
workflow at least 5000. Raise concurrency toward 50 after checking progress
and quota. The checked-in `google_batch.conf` is a cold, bounded benchmark
configuration and must not replace that server configuration. Client options
do not reconfigure the server. Use the existing locality guard with the actual
server execution root, Batch region, inputs and cohort options.

Run two full-depth libraries first as part of the intended 100, including one
typical library and one larger prior RSEM case. Check gene/isoform outputs,
expression metadata and [task profiles](task-profiling.md), then expand across
read depth, prior BAM/runtime extremes and batches. Reuse the first two through
verified cache hits. Preserve inputs/options and source revision. An independent
sample can finish after another fails; a full merge still requires repair of
the missing sample. Capture failed attempts as well as successes before cleanup.

Keep the N1 family and the eligible predefined RNA-QC upgrade for this calibration. A price-only comparison found
E2 predefined Spot shapes worth a separate small matched test; changing family
at the same time would confound RAM calibration. Predefined shape selection
overrides CPU/RAM runtime requests, so never apply one fixed machine type to
the entire workflow or bypass dynamic RSEM growth.
