# Rat GRCr8 / Ensembl 116

Use `-a rat -v rn8_v116` with
`--runtime-profile config/backends/gcp/runtime-rat-v116-full-candidate-v1.json`
and `--star-disk-type SSD`. The release profile pins the same tool/container
versions as human GENCODE v50, including STAR 2.7.11b, RSEM 1.3.3, Bowtie2 2.5.5,
Subread 2.1.1, Samtools 1.24 and Picard 3.5.0. `rn8` continues to mean Ensembl 115.
Hand-edited inputs must also set `rnaseq_pipeline.reference_release` to
`rn8_v116`; this controls the reference label and rat RSEM memory sizing.
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
Floors are 40 GiB STAR RAM / 120 GB SSD, 11 GiB UMI RAM / 30 GB SSD,
18 GiB RSEM RAM / 30 GB scratch, and 3 GiB RNA-QC RAM. Across 60 full-depth
libraries, including five tissues and 54 liver samples, peak working RAM was
33.7/9.4/15.1/1.9 GiB respectively. STAR, UMI and RSEM allocations provide
17–20% above those maxima; RNA QC rounds up to whole GiB. MarkDuplicates retains
36 GiB against a 31.0-GiB peak. Working RAM excludes reclaimable inactive file
cache; reduced cache space may affect runtime.

Rat RSEM RAM is `max(configured floor, ceil(0.5 + 4 × input BAM GiB))`.
The same rule applies to molecule and optional all-read expression. Other
releases retain their existing formula. Shared scratch growth from post-trim
pairs and the BAM actually entering each step remains active. Calibration covers
27.2–41.4 million post-trim pairs and 0.63–4.31 GiB molecule transcriptome BAMs;
larger inputs and optional all-read fits remain extrapolations. STAR and UMI RAM
are configured allocations, not read-depth growth formulas.

These reduced allocations require validation on subsequent production samples.
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
