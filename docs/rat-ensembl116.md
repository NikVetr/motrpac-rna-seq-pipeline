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

## Resource calibration and validation

The candidate requests E2 for STAR, UMI, RSEM, RNA QC and the benchmarked
preprocessing/counting/QC tasks listed in [the E2 policy](v50-cohort-calibration.md).
Floors are
48 GiB STAR RAM / 120 GB SSD, 24 GiB UMI RAM / 80 GB SSD, 24 GiB RSEM RAM /
60 GB scratch, and 8 GiB RNA-QC RAM. Shared growth from post-trim pair counts and
the BAM actually entering each step remains active. These are buffered starting
allocations, not measured rat minima. Six full-depth libraries completed all
80 calls on their first attempts, with 28.5–41.4 million post-trim pairs across
five tissues. Peak working RAM was 33.7 GiB for STAR, 9.4 GiB for UMI, 15.1 GiB
for RSEM and 1.9 GiB for RNA QC. Peak scratch was 89.9/21.8/18.9/3.5 GiB,
respectively. Keep these tested allocations for broader rat calibration;
the UMI and RNA-QC floors offer the clearest subsequent RAM reductions.

All 88 local unit tests pass, along with MiniWDL, WOMtool input validation and
rendered WDL resource/expression checks. These cover contig conversion,
mapped-QC aliases, matched release inputs and shared scientific/output contracts.
The two-library, 100,000-pair cloud gate passed all 30 calls, including an explicit
missing-I1 case and optional all-read outputs. Gene/transcript matrix values match
the raw RSEM files exactly; expression metadata records the expected UMI policy.
Independent forward/reverse featureCounts checks support forward strandedness in
both muscle and blood. Full-depth PASS1B libraries spanning gastrocnemius,
blood, cortex, liver and white adipose also pass exact raw-to-matrix and UMI
metadata checks. Their worker cost estimate is $5.06 total, or $0.84/library on
demand ($0.52/library repriced at uninterrupted Spot rates).

An eleven-call E2 comparison reused one full-depth muscle library's inputs and
the six-library merge inputs. All 32 scientific output checks passed. Seven
preprocessing/counting/QC calls cost 24–37% less individually on E2; chromosome
and QC reporting plus the two merges cost 18–40% more, adding $0.007 combined.
The E2 policy therefore selects the seven winners, retaining backend selection
for those four calls. These results are one matched replay per task, not a
hardware guarantee or a human-cohort benchmark. Scientific commands are unchanged.

Evidence is archived under
`gs://motrpac-rnaseq-modernization-us-west1/rat-ensembl116-pilot-20260912/controller-final/`:
full workflow `9dc8431f-796d-455c-8961-16f513537be8`, E2 comparison
`3211459d-a378-448c-bd17-e6ea74f27101`. Cost estimates use captured September 12,
2026 us-west1 rates and Batch startup/run durations, including boot/scratch;
they exclude controller, teardown, retained storage and external downloads.
All pilot workers were removed and the controller stopped. Rat compatibility
MultiQC and warm-cache recovery were not separately cloud-tested.
