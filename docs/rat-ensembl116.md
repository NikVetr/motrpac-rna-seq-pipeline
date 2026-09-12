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
guard validates execution placement. Production options can retain caching and
`ContinueWhilePossible`; a cold benchmark explicitly disables caching and retries.

## Initial resource calibration

The candidate requests E2 for STAR, UMI, RSEM and RNA QC. Initial floors are
48 GiB STAR RAM / 120 GB SSD, 24 GiB UMI RAM / 80 GB SSD, 24 GiB RSEM RAM /
60 GB scratch, and 8 GiB RNA-QC RAM. Shared growth from post-trim pair counts and
the BAM actually entering each step remains active. These are buffered starting
allocations, not measured rat minima.

All 88 local unit tests pass, along with MiniWDL, WOMtool input validation and
rendered WDL resource/expression checks. These cover contig conversion,
mapped-QC aliases, matched release inputs and shared scientific/output contracts.
The two-library, 100,000-pair cloud gate passed all 30 calls, including an explicit
missing-I1 case and optional all-read outputs. Gene/transcript matrix values match
the raw RSEM files exactly; expression metadata records the expected UMI policy.
Independent forward/reverse featureCounts checks support forward strandedness in
both muscle and blood. Six full-depth PASS1B libraries spanning gastrocnemius,
blood, cortex, liver and white adipose are under calibration. An eleven-call E2
comparison reuses retained inputs for the smaller steps. Full-depth performance
and resource calibration are pending.
