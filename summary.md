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
enabled by default when matched I1 reads are present. Conventional all-read
matrices are always emitted. The generator exposes an explicit legacy switch
that omits molecule matrices without changing the submission interface.

FASTQ QC, contamination screening, alignment QC, and UMI QC are independently
selectable. Native tool reports feed the stable QC table directly. The three
contamination screens can share one worker and one deterministic post-trim read
sample, while full-depth screening remains available. An opt-in compatibility
mode publishes the legacy pre- and post-alignment MultiQC archives without
changing the QC table or expression outputs. It uses the immutable production
MultiQC 1.6 image and requires both FastQC groups and alignment QC.

Complete runtime profiles provide explicit CPU, memory, disk size, and STAR
disk-class settings for bounded GCP canaries. The GCP support layer includes a
concurrency guard, read-only preflight, pinned Cromwell Batch configuration,
resource monitoring, immutable evidence capture, and attempt-aware cost
summarization. Generated evidence and rendered benchmark reports are analysis
artifacts and are not part of the production repository.

The focused 60-test suite covers input validation, release/runtime profiles,
WDL I/O contracts, native QC parsing, contamination sampling, directional UMI
grouping, molecule-expression construction, and the GCP monitoring/cost
contracts. The production execution tree also passes MiniWDL and WOMtool 91
validation under OpenJDK 21.

End-to-end acceptance covers the intended human-v47 default graph with 17 calls
and the opt-in MultiQC compatibility graph with 19 calls. The retained v39 and
rat configurations, no-I1 policy, QC switches, and legacy all-read-only mode
have focused contract coverage, but every cross-product of those alternate
settings has not been run as a separate integration workflow. The exact
high-input resource vector remains outside end-to-end acceptance.

The default human-v47 graph passes complete 100k-pair integration canaries both
locally and through Cromwell 92/GCP Batch. The exact on-demand cloud execution
completed all 17 calls on their first attempt, emitted all 12 expected
top-level outputs, stayed within the controller-plus-three-worker cap, and had
a modeled worker cost of $0.15889. All eight conventional and molecule-level
matrices contain 78,932 unique, finite, nonnegative genes and are byte-identical
between local and cloud. The QC CSV, UMI metrics, and contamination manifest
are also byte-identical; only a temporary compressed-BAM byte-size provenance
field differs. Cromwell's Cloud SDK helper is pinned to an immutable official
Python-equipped image. Local MiniWDL ignores WDL disk requests, and the exact
150-GiB high-input STAR profile has not yet processed a full large sample. The
dedicated cloud controller is stopped and no Batch worker remains running.

The MultiQC-enabled 100k canary completed all 19 calls locally and through
Cromwell 92/GCP Batch. Cloud workflow
`841a7f5a-7780-43a8-9cba-0b518e6a5a79` completed every call on its first
attempt, emitted all 14 top-level outputs, and stayed within the controller plus
three-worker cap. Both report archives contain a nonempty HTML report and the
expected module data: pre-alignment reporting detected Cutadapt and four
FastQC reports; post-alignment reporting detected RSEM, featureCounts, STAR,
Cutadapt, two FastQC reports, Picard MarkDuplicates, and Picard RNA metrics.
All eight expression matrices and the QC CSV are byte-identical to the accepted
cloud canary without MultiQC. The controller is stopped and no Batch worker
remains running. Registration and output ingestion through the external
operator GUI/Caper deployment remain a separate handoff check.
