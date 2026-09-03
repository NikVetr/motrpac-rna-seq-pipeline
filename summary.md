# Project summary

This repository implements the MoTrPAC paired-end RNA-seq workflow in WDL for
Cromwell/Caper. The main scatter workflow validates paired FASTQ and optional
I1 inputs, attaches UMIs, trims adapters, aligns with STAR, quantifies with
RSEM and featureCounts, runs selectable QC branches, and gathers cohort-ready
matrices and QC outputs.

Human runs support the historical GENCODE v39 configuration and an immutable
GENCODE v47 release profile whose references, index-builder/runtime versions,
and container digests are validated as one unit. Rat rn6, rn7, and rn8 inputs
retain their existing configurations.

Directional UMI grouping and molecule-level RSEM and featureCounts matrices are
enabled by default when matched I1 reads are present. These matrices use the
canonical output names. Conventional all-read quantification is skipped unless
the operator requests secondary `all_read_*` matrices; an explicit legacy
switch instead makes historical all-read expression canonical.

FASTQ QC, contamination screening, alignment QC, and UMI QC are independently
selectable. Native tool reports feed the stable QC table directly. The three
contamination screens can share one worker and one deterministic post-trim read
sample, while full-depth screening remains available. An opt-in compatibility
mode publishes the legacy pre- and post-alignment MultiQC archives without
changing the QC table or expression outputs. It uses the immutable production
MultiQC 1.6 image and requires both FastQC groups and alignment QC.

Complete runtime profiles provide explicit CPU, memory, disk floors, and STAR
disk-class settings for bounded GCP canaries. Cutadapt exposes the exact
surviving read-pair count, from which the workflow selects a buffered 90-,
120-, 150-, 180-, 200-, 250-, 300-, or 400-GB STAR scratch tier independently
for every sample; the selected profile value remains a minimum. A single input
JSON therefore supports heterogeneous sample sizes. Exact include/exclude
sample manifests support bounded pilots followed by nonoverlapping cohort
runs. The GCP support layer includes a concurrency guard, read-only preflight,
pinned Cromwell Batch configuration, resource monitoring, immutable evidence
capture, and attempt-aware cost summarization. With no `cpuPlatform` override,
Cromwell provisions N1 custom workers; cost summaries therefore default to a
frozen N1 manifest and require an explicit family-matched manifest for N2
evidence. Generated evidence and rendered benchmark reports are analysis
artifacts and are not part of the production repository.

The focused 64-test suite covers input validation, release/runtime profiles,
WDL I/O contracts, native QC parsing, contamination sampling, directional UMI
grouping, molecule-expression construction, and the GCP monitoring/cost
contracts. The production execution tree also passes MiniWDL and WOMtool 91
validation under OpenJDK 21.

The current human-v47 graph, retained v39 and rat configurations, no-I1 policy,
QC switches, default molecule-expression policy, optional all-read branch, and
legacy all-read-only mode have focused contract coverage, but every
cross-product has not been run as a separate integration workflow. The current
per-sample disk-tier and canonical-output policies require one bounded
multi-sample production-interface pilot before cohort release.

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
