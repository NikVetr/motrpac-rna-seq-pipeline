# Project summary

This repository implements the MoTrPAC paired-end RNA-seq workflow in WDL for
Cromwell/Caper. The main scatter workflow validates paired FASTQ and optional
I1 inputs, attaches UMIs, trims adapters, aligns with STAR, quantifies with
RSEM and featureCounts, runs selectable QC branches, and gathers cohort-ready
matrices and QC outputs.

Human runs support the historical GENCODE v39 configuration and immutable
GENCODE v47 and v50 release profiles whose references, index-builder/runtime versions,
and container digests are validated as one unit. Rat rn6, rn7, and rn8 inputs
retain their existing reference configurations.
The v50 primary-assembly references contain 646,577 transcripts, use the same
STAR 2.7.11b and RSEM 1.3.3 tools as v47, and reside in private us-west2 storage.
Published v47 resource mappings remain specific to their measured annotation.

Directional UMI grouping and molecule-level RSEM and featureCounts matrices are
enabled by default when matched I1 reads are present. These matrices use the
canonical output names. Conventional all-read quantification is skipped unless
the operator requests secondary `all_read_*` matrices. The explicit
`--all-read-expression-only` switch instead makes all-read expression canonical,
without changing strandedness; its former spelling remains an alias. Combining
it with `--no-index` omits I1 inputs and UMI QC. Opt-in `--allow-missing-umis`
uses molecule expression for samples with I1 and all-read expression for samples
without it. Strict mode still rejects missing UMIs. Canonical expression metadata
records sample order, reference release, UMI availability, `not_deduplicated`,
expression mode and skip reason; it can be joined to study/QC metadata.

Raw RSEM gene and isoform results are canonical outputs, preserving effective
lengths and transcript-to-gene mappings. RSEM produces both in the same task.
Gene and transcript count/TPM/FPKM matrices are separate outputs. The transcript
merge streams matching rows across samples and rejects inconsistent transcript
order or gene mappings, avoiding cohort-sized in-memory transcript matrices.

FASTQ QC, contamination screening, alignment QC, and UMI QC are independently
selectable. Native tool reports feed the stable QC table directly. The three
contamination screens can share one worker and one deterministic post-trim read
sample, while full-depth screening remains available. An opt-in compatibility
mode publishes the legacy pre- and post-alignment MultiQC archives without
changing the QC table or expression outputs. It uses the immutable production
MultiQC 1.6 image and requires both FastQC groups and alignment QC.

Complete runtime profiles provide explicit CPU, memory, disk floors, and
operator-selected STAR and UMI disk classes for bounded GCP canaries. Cutadapt
exposes the exact surviving read-pair count, from which the workflow selects a
buffered 90-, 120-, 150-, 180-, 200-, 250-, 300-, or 400-GB STAR scratch tier
independently for every sample. UMI scratch is independently raised to twice
the combined STAR BAM size plus 15 GB. The selected profile values remain
minimums, including the benchmarked 80-GB full-depth UMI floor. UMI scratch
defaults to SSD with an explicit HDD override; STAR remains independently
selectable and defaults to HDD. A single input
JSON therefore supports heterogeneous sample sizes. Exact include/exclude
sample manifests support bounded pilots followed by nonoverlapping cohort
runs. Both cohort merges request at least three times their total input GiB
plus 10 GiB, rounded up, preserving the configured disk request as a floor.
Merge inputs are symbolic links to localized files; the validated disk margin
is retained. Merge RAM grows by 4 GB per 75 libraries, rounded up, with the
configured value as a floor. Molecule RSEM RAM and scratch also grow from
transcriptome BAM bytes, preserving configured floors. CPU counts are explicit.
The merged QC table retains the legacy
pipeline-derived covariate contract; participant, visit, treatment, demographic,
batch, and RIN metadata are joined from study records by sample ID using
`scripts/prepare_sample_metadata.py`, preserving expression-column order and
unscaled values. A launch guard requires regional execution storage matching
the Batch region and worker zones, while reporting remote source buckets.
An explicit controller-stop helper requires complete local evidence first.
The GCP support layer includes a concurrency guard, read-only preflight,
pinned Cromwell Batch configuration, resource monitoring, immutable evidence
capture, and attempt-aware cost summarization. With no `cpuPlatform` override,
Cromwell provisions N1 custom workers; cost summaries therefore default to a
frozen N1 manifest and require an explicit family-matched manifest for N2
evidence. Generated evidence and rendered benchmark reports are analysis
artifacts and are not part of the production repository.

The focused 74-test suite covers input validation, release/runtime profiles,
WDL I/O contracts, native QC parsing, contamination sampling, directional UMI
grouping, molecule-expression construction, and the GCP monitoring/cost
contracts. The production execution tree also passes MiniWDL and WOMtool 91
validation under OpenJDK 21.
`tests/check_wdl_resources.py` additionally evaluates resource expressions and
checks both rendered merge commands against byte-identical fixture outputs.
The cohort resource review covers all 4,455 successful per-sample calls from
297 libraries; its full RSEM memory peak is 38.26 GiB, supporting retention of
the 40-GB floor. A three-call GCP canary verified RSEM growth to 20-GB RAM and
15-GB scratch and both two-library merges, with byte-identical outputs, one
worker at a time, and no retries. The cohort/canary evidence does not establish
large-sample speedups or benefits from smaller SSDs.

The current human-v47 graph, retained v39 and rat configurations, no-I1 policy,
QC switches, default molecule-expression policy, optional all-read branch, and
legacy all-read-only mode have focused contract coverage, but every
cross-product has not been run as a separate integration workflow. The human
v47 cohort completed all 4,440 per-sample calls. A combined 297-library gather,
including the separately completed missing sample, succeeded with 68 GB of
dynamically sized scratch and an explicit 16-GB RAM request. All four matrices
passed gene/sample and numeric checks; merged QC matched every source row.

The current canonical-output graph completed a 100k-pair local canary in 16
calls with MultiQC enabled. It did not schedule conventional all-read
featureCounts, RSEM, or gather calls; every optional `all_read_*` output was
empty. The canonical featureCounts and RSEM count matrices were byte-identical
to the directional-UMI matrices from the preceding accepted local canary. The
Cutadapt task reported exactly 99,979 surviving pairs to the STAR sizing
expression. Local MiniWDL does not provision WDL disks; cloud runtime evidence
is required to verify the selected disk tiers.

The preceding dual-expression human-v47 graph passes complete 100k-pair
integration canaries both locally and through Cromwell 92/GCP Batch. The exact
on-demand cloud execution completed all 17 calls on their first attempt,
emitted all 12 expected top-level outputs, stayed within the
controller-plus-three-worker cap, and had a modeled worker cost of $0.15889.
All eight conventional and molecule-level matrices contain 78,932 unique,
finite, nonnegative genes and are byte-identical between local and cloud. The
QC CSV, UMI metrics, and contamination manifest are also byte-identical; only a
temporary compressed-BAM byte-size provenance field differs. Cromwell's Cloud
SDK helper is pinned to an immutable official Python-equipped image. Local
MiniWDL ignores WDL disk requests. Exact revision
`0a89dd15b1be05b27781902d353618003e639f7f` processed sample `11076050401`
(48,078,786 read pairs) with the 150-GiB high-input STAR profile, forward
strandedness, directional molecule expression, sampled contamination QC, and
MultiQC enabled. Workflow `9f950a98-17fa-4278-84de-c295c43cfd81` completed all
19 calls on their first on-demand attempt and emitted all 14 expected top-level
outputs. The operator interface also completed a full-depth dual-expression
human-v47 run.

The MultiQC-enabled 100k canary completed all 19 calls locally and through
Cromwell 92/GCP Batch. Cloud workflow
`841a7f5a-7780-43a8-9cba-0b518e6a5a79` completed every call on its first
attempt, emitted all 14 top-level outputs, and stayed within the controller
plus three-worker cap. Both report archives contain a nonempty HTML report and
the expected module data: pre-alignment reporting detected Cutadapt and four
FastQC reports; post-alignment reporting detected RSEM, featureCounts, STAR,
Cutadapt, two FastQC reports, Picard MarkDuplicates, and Picard RNA metrics.
All eight expression matrices and the QC CSV are byte-identical to the accepted
cloud canary without MultiQC. No benchmark controller or Batch worker is
intended to remain running after evidence capture.
