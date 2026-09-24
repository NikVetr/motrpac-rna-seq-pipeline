# Rat GRCr8 / Ensembl 116

Use `-a rat -v rn8_v116` with
`--runtime-profile config/backends/gcp/runtime-rat-v116-full-candidate-v1.json`
and `--star-disk-type SSD`. The release profile pins the same tool/container
versions as human GENCODE v50, including STAR 2.7.11b, RSEM 1.3.3, Bowtie2 2.5.5,
Subread 2.1.1, Samtools 1.24 and Picard 3.5.0. `rn8` continues to mean Ensembl 115.
Hand-edited inputs must also set `rnaseq_pipeline.reference_release` to
`rn8_v116`; this controls the reference label and rat UMI/RSEM memory sizing.
Pulling the repository does not update already-generated input JSONs.

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
Floors are 38 GiB STAR RAM / 120 GB SSD, 8 GiB UMI RAM / 30 GB SSD,
6 GiB molecule-RSEM RAM / 30 GB scratch, and 3 GiB RNA-QC RAM. STAR's allocation
retains 12.8% headroom over the largest measured working set; MarkDuplicates
retains 36 GiB against a 32.5-GiB peak. Working RAM excludes reclaimable inactive
file cache; reduced cache space may affect runtime.

UMI RAM grows as `ceil(8.1 + 1.2 × genomic BAM GiB)` and molecule-RSEM RAM as
`ceil(2 + 4 × input BAM GiB)`, preserving larger configured floors. The rules
use OLS predictions plus the larger of the held-out 99th-percentile residual
and three residual standard deviations, then 15% headroom and upward rounding.
Validation groups animals and holds out whole tissues. The calibration includes
1,506 UMI and 1,505 RSEM calls across 19 tissues, with a separate deeper panel
covering up to 106 million post-trim pairs. RSEM memory is strongly predicted
by input BAM size; UMI's larger residual allowance accommodates tissue variation.
Failed-call peaks are lower bounds, so coverage of completed calls does not
guarantee recovery of a memory-exhausted sample or 99% coverage of new tissues.

UMI requests two CPUs; E2 selection raises this when RAM exceeds 16 GiB.
RSEM retains ten threads and an 18-GiB minimum for optional all-read expression.
Shared scratch growth and disk floors remain active to preserve I/O throughput.
The smaller VM shapes require runtime/cost confirmation on subsequent samples.
Generate new inputs with the updated runtime profile and retain monitoring;
do not resubmit completed cohorts solely to adopt provisioning changes. Changed
task inputs, including the RSEM release selector, can invalidate call-cache hits.

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
