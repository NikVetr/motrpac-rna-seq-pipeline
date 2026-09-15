# Project summary

The MoTrPAC paired-end RNA-seq pipeline runs through WDL and Cromwell/Caper.
`wdl/rnaseq_pipeline_scatter.wdl` scatters over samples, attaches optional I1
UMIs, trims adapters, aligns with STAR, quantifies with RSEM and featureCounts,
collects QC, and merges expression matrices and sample metadata.

## References and expression

Human GENCODE v47/v50 and rat GRCr8/Ensembl 116 use matched reference and image
profiles in `config/release-profiles/`. Reference manifests record source
checksums, index-building parameters and object generations; container images
are pinned by digest. The rn6, rn7, rn8/Ensembl 115 and human v39 configurations
remain selectable. Reference selection is independent of execution region.

Counting is forward-stranded. Directional UMI molecule expression is canonical
when I1 reads are present. Missing UMIs fail unless `--allow-missing-umis` is
enabled; permitted missing-I1 samples use all-read expression and record
`not_deduplicated=1` and a skip reason in expression metadata.
`--all-read-expression-only` selects all-read counting for every sample;
`--retain-all-read-expression` adds a secondary all-read pass to molecule
counting. Neither option changes strandedness.

RSEM emits gene and isoform results in the same calculation. Canonical outputs
retain both raw result types, effective lengths and transcript-to-gene mappings,
plus gene/transcript count, TPM and FPKM matrices. Transcript merging streams
rows and checks transcript order and gene mappings. Study covariates are joined
by sample ID using `scripts/prepare_sample_metadata.py`.

## QC and resources

Native tool reports feed the QC table. FastQC, contamination, alignment and UMI
QC are selectable; optional MultiQC archives require both FastQC groups and
alignment QC. Combined contamination screening shares one worker and an optional
deterministic post-trim read sample across globin, rRNA and PhiX screens.

Runtime profiles set resource floors. STAR scratch grows with post-trim pairs,
UMI scratch with combined STAR BAM bytes, and RSEM RAM/scratch with the
transcriptome BAM entering that task. Merge resources grow with library count
and input bytes. Human v50 applies a separate STAR scratch multiplier.
The human v50 and rat v116 profiles select E2 for supported processing tasks;
one-CPU reporting and merges retain backend selection. See
[cohort provisioning](docs/cohort-provisioning.md) for rules and deployment.

## Operation and validation

`scripts/make_json_rnaseq.py` validates paired inputs, release/runtime profiles,
and exact include/exclude sample lists. GCP helpers check locality, monitor
resources, capture attempt-level evidence and estimate costs from frozen rates.
Execution storage and workers must share a region; source references and FASTQs
require explicit read permissions and placement decisions. Cohort submissions
use persistent call caching and `ContinueWhilePossible`; the separate canary
backend limits concurrency and disables caching for controlled measurements.

Tests cover input and reference contracts, QC parsing, UMI grouping and molecule
propagation, gene/isoform merging, profiling and failure accounting. Run
`python3 -m unittest discover -s tests`; in a MiniWDL environment, also run
`tests/check_wdl_resources.py` and `tests/check_wdl_expression_policy.py`.
Validate the workflow before deployment. GCP task-specific machine selection
requires Cromwell 92. App deployments must synchronize their bundled workflow,
input settings and backend with the selected pipeline release.
