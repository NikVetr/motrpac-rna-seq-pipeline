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
sample, while full-depth screening remains available.

Complete runtime profiles provide explicit CPU, memory, disk size, and STAR
disk-class settings for bounded GCP canaries. The GCP support layer includes a
concurrency guard, read-only preflight, pinned Cromwell Batch configuration,
resource monitoring, immutable evidence capture, and attempt-aware cost
summarization. Generated evidence and rendered benchmark reports are analysis
artifacts and are not part of the production repository.

Tests under `tests/` cover input validation, release/runtime profiles, WDL I/O
contracts, native QC parsing, contamination sampling, directional UMI grouping,
molecule-expression construction, and the GCP monitoring/cost contracts.
