# Rat GRCr8 / Ensembl 116

Use `-a rat -v rn8_v116` with
`--runtime-profile config/backends/gcp/runtime-rat-v116-full-candidate-v1.json`
and `--star-disk-type SSD`. The release profile pins the same tool/container
versions as human GENCODE v50, including STAR 2.7.11b, RSEM 1.3.3, Bowtie2 2.5.5,
Subread 2.1.1, Samtools 1.24 and Picard 3.5.0. `rn8` continues to mean Ensembl 115.

The shared workflow performs directional UMI deduplication and forward-stranded
RSEM/featureCounts quantification by default. Canonical outputs include raw RSEM
gene and isoform results, gene/transcript matrices and expression metadata.
`--allow-missing-umis` explicitly permits mixed I1 availability and records the
skip reason per sample; `--retain-all-read-expression` adds optional all-read
outputs. Rat uses the same combined, sampled contamination QC, native QC parsers,
cohort merge sizing, step-input profiling and E2 selection as human runs.

## Reference contract

[`rat-grcr8-ensembl-v116.json`](../config/references/rat-grcr8-ensembl-v116.json)
records Ensembl download checksums, normalized FASTA/GTF hashes, build versions,
payload hashes and immutable GCS object generations. The source is the complete
[Ensembl 116 GRCr8 annotation](https://jun2026.archive.ensembl.org/Rattus_norvegicus/Info/Annotation)
and unmasked toplevel assembly GCA_036323735.1: 2,849,597,492 bases, 43,360 genes
and 95,472 transcripts. Numeric/X/Y chromosomes receive `chr` prefixes, MT becomes
chrM, and scaffold names remain unchanged. Mapped-QC accepts both naming styles.

STAR uses the normalized FASTA/GTF with `sjdbOverhang=100`; RSEM uses exactly the
same inputs. UCSC gtfToGenePred kent-v479 generates an 11-column refFlat with one
row per transcript. The rat globin/rRNA and shared PhiX assay target sequences
are retained and their indexes rebuilt with Bowtie2 2.5.5; these curated QC
targets are not derived from the new Ensembl annotation.

References are checksum-addressed private objects in
`gs://motrpac-rnaseq-modernization-us-west1`. Workers need read access to this
bucket and the raw FASTQs. Use us-west1 execution storage and worker zones for
the PASS1B rat inputs, which are already in us-west1. The submission locality
guard validates execution placement. Use caching and `ContinueWhilePossible`
for cohorts; configure cache and retry settings explicitly for cold benchmarks.

## Resource calibration and validation

The candidate requests E2 for STAR, UMI, RSEM, RNA QC and the benchmarked
preprocessing/counting/QC tasks listed in [the E2 policy](v50-cohort-calibration.md).
Floors are
48 GiB STAR RAM / 120 GB SSD, 24 GiB UMI RAM / 80 GB SSD, 24 GiB RSEM RAM /
60 GB scratch, and 8 GiB RNA-QC RAM. Shared growth from post-trim pair counts and
the BAM actually entering each step remains active. These are buffered starting
allocations, not measured rat minima. A six-library, five-tissue calibration
spanned 28.5–41.4 million post-trim pairs. Peak working RAM was 33.7/9.4/15.1/1.9
GiB for STAR/UMI/RSEM/RNA-QC; peak scratch was 89.9/21.8/18.9/3.5 GiB.
Retain these buffered floors while collecting broader cohort measurements.

Validation covers contig conversion, mapped-QC aliases, matched references,
forward strandedness, mixed I1 availability, optional all-read outputs, and exact
raw-RSEM-to-matrix agreement. Full-depth coverage includes muscle, blood,
cortex, liver and white adipose. E2 task comparisons support the shared selection
policy; performance remains workload-dependent. Rat compatibility MultiQC and
warm-cache recovery require deployment validation.

Calibration evidence is retained under
`gs://motrpac-rnaseq-modernization-us-west1/rat-ensembl116-pilot-20260912/controller-final/`.
Follow [cohort deployment](cohort-provisioning.md#deployment-and-locality) for
the two-sample acceptance gate, broader sampling, profiling and output retention.
