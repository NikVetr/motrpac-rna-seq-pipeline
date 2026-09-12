# Initial human v50 cohort calibration

Use `config/backends/gcp/runtime-human-v50-full-candidate-v1.json` for the
initial approximately 100 libraries. It preserves the v47 full-depth profile's
threads and disk floors, changes UMI/RSEM/RNA-QC RAM floors to 36/32/12 GB,
and keeps STAR at 12 threads/72 GB. Shared BAM-size growth remains active.
It enables `rnaseq_pipeline.use_e2` for STAR, UMI, both RSEM expression modes
and RNA-QC. Other tasks retain their existing backend sizing.
The workflow multiplies STAR's post-trim read-pair scratch tiers by 1.30 only
when `reference_release` is `gencode_v50`, rounding up and respecting larger
explicit floors. v47 and rat resource rules are unchanged.

E2 custom machines retain the requested RAM and use the smallest even vCPU
count satisfying both tool threads and the 8-GiB/vCPU memory limit. RSEM uses
its effective BAM-scaled RAM, including when growth exceeds the configured
floor. For the tested libraries, STAR uses 12/72 CPU/GiB, UMI 6/36, RSEM
10/32–56 and RNA-QC 2/12. Tool threads and disk formulas are unchanged.
Do not apply one fixed machine type to the entire workflow or clamp growing
RAM requests. Requests beyond E2's supported sizes fail allocation and need an
explicit profile/family decision. The locality guard rejects simultaneous
predefined N1 selection and CPU-platform/machine overrides. Use Cromwell 92,
which supports these conditional `gcp` runtime objects. For the N1 comparison
mode, disable `use_e2` and optionally enable `prefer_predefined_n1` in us-west2.

These are buffered pilot candidates, not validated resource minima. Four
completed v50 RSEM calls used 7.38/27.82/31.85/35.71 GiB working memory and
would request 32/48/52/56 GB. Incomplete UMI tails reached 30.94 GiB working
memory; completed RNA-QC reached approximately 10 GiB. Larger v50 tails and
the cohort transcript merge still need measurements at these allocations.
All 12 heavy-task E2 calls completed across three full-depth v50 libraries.
Their RSEM working-memory peaks were 7.38/30.39/34.51 GiB; STAR was below
33 GiB, UMI below 16 GiB and RNA-QC below 9.5 GiB. Every allocation retained
disk headroom. RSEM gene and isoform tables were byte-identical to matched N2
outputs; available UMI counters, STAR metrics/junctions and Picard metrics
also agreed. This supports a staged 100-library calibration, not further
resource reductions or a claim that larger tails have already been tested.

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
expression metadata and [task profiles](task-profiling.md), then expand to ten
including high-depth/resource extremes, and finally all 100. Reuse completed
libraries through verified cache hits. Preserve inputs/options and source revision. An independent
sample can finish after another fails; a full merge still requires repair of
the missing sample. Capture failed attempts as well as successes before cleanup.

Keep the candidate RAM/disk relationships fixed during this E2 calibration.
Select approximately 80 representative libraries stratified by batch, visit
and depth, plus 20 anchors/resource extremes from the original 297. Retain
the selection reasons and historical metrics separately from newly generated
QC. Reweight the deliberately enriched panel before estimating cohort-average
costs. No additional all-read RSEM pass is needed to produce isoform results.
