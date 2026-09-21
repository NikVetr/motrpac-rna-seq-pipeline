# Human v50 cohort calibration

Use `config/backends/gcp/runtime-human-v50-full-candidate-v1.json` for the
cohort. STAR uses 12 threads/56 GB RAM; UMI, RSEM and RNA-QC RAM floors
are 20, 32 and 12 GB. STAR/UMI/RSEM scratch floors are 120/80/60 GB.
[Shared resource rules](cohort-provisioning.md) scale v50 STAR scratch with
post-trim pairs and UMI/molecule RSEM memory with task input BAM sizes.
v47, rat and all-read RSEM memory rules are unchanged.

## E2 policy and calibration scope

`rnaseq_pipeline.use_e2` selects E2 for STAR, UMI, both RSEM expression modes,
RNA-QC, FastQC, UMI attachment, Cutadapt, featureCounts, combined contamination
QC and MarkDuplicates. One-CPU chromosome/QC reporting, merges and optional
legacy QC retain backend sizing.

E2 custom machines retain requested RAM and use the smallest even vCPU count
satisfying both tool threads and the 8-GiB/vCPU memory limit. RSEM uses its
effective BAM-scaled RAM. Tool threads need not equal allocated vCPUs.
Do not apply a fixed machine type or CPU platform to the entire workflow.
Requests beyond E2's supported sizes require an explicit profile/family decision.
Use Cromwell 92. For N1 comparisons, disable `use_e2`; the optional
`prefer_predefined_n1` policy is restricted to us-west2.

These buffered allocations apply to full-depth human muscle libraries; they
are not universal resource bounds for other tissues or larger inputs. Higher
configured floors remain effective. Existing input JSONs retain their saved
floors; regenerate inputs with this runtime profile to use its current settings.

## Inputs and submission

Choose approximately 100 intended production libraries spanning batch, visit,
depth and known UMI/RSEM resource extremes. Record selection reasons and reweight
an extreme-enriched panel when estimating whole-cohort costs. Generate inputs
from each FASTQ directory using an exact sample list:

```bash
python3 scripts/make_json_rnaseq.py \
  -g gs://BUCKET/BATCH/fastq_raw -o input_json -r v50_calibration \
  -a human -v gencode_v50 -n 1 -p PROJECT \
  --sample-list selected_samples.txt \
  --runtime-profile config/backends/gcp/runtime-human-v50-full-candidate-v1.json \
  --star-disk-type SSD --umi-dup-disk-type SSD \
  --combine-contamination-qc --contamination-qc-pairs 1000000 -i
```

Add `--allow-missing-umis` for mixed I1 availability; those samples are explicitly
marked as not deduplicated. Default counting produces molecule-level gene and
isoform results in one RSEM pass. A secondary all-read pass is optional.

Copy `config/backends/gcp/workflow-options-cohort.example.json` and set zones
and monitoring-script placement for the actual execution region. It enables
cache reads/writes, `ContinueWhilePossible`, and `maxRetries: 0` for command
failures. The workflow defaults to one Spot attempt before on-demand fallback;
an explicit `rnaseq_pipeline.num_preemptible_attempts` overrides this
(`0` for on-demand only, `2` for two Spot attempts). Command retries and
preemption retries are separate; neither raises RAM after an OOM.

Follow the server and app requirements in [cohort deployment](cohort-provisioning.md#deployment-and-locality).
Run two full-depth libraries from the intended cohort first, including a typical
sample and a resource extreme. Verify gene/isoform outputs, expression metadata
and [task profiles](task-profiling.md), then expand to the remaining cohort.
Preserve accepted results, source revision and input/options manifests.

Keep candidate allocations fixed within the calibration batch. Include failed
attempts and retries in the resource/cost analysis, and repair failed samples
before the full gather. Provisioning changes can follow in a separate release
without regenerating accepted results.
