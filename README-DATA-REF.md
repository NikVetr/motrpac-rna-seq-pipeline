# Preparing Reference Files for the RNA-seq Pipeline

This document describes how to prepare reference genome files for running the MoTrPAC RNA-seq pipeline with a new organism or genome build.

The example below uses the rat rn8 (GRCr8) assembly with Ensembl release 115 annotations. The same process was used for rn6 (Ensembl 96) and rn7 (Ensembl 108).

---

## Overview

The RNA-seq pipeline requires **7 reference files**:

| # | File | Description | How to Create |
|---|------|-------------|---------------|
| 1 | **Genome FASTA** | Reference genome sequence | Download from Ensembl FTP, add `chr` prefix |
| 2 | **GTF file** | Gene annotation in GTF format | Download from Ensembl FTP, add `chr` prefix |
| 3 | **STAR index** | Genome index for STAR aligner | Build with WDL workflow |
| 4 | **RSEM reference** | Reference for RSEM quantification | Build with WDL workflow |
| 5 | **refFlat file** | Gene intervals for Picard QC | Generate from chr-prefixed GTF |
| 6 | **Globin index** | Bowtie2 index for globin contamination | Build or reuse existing |
| 7 | **rRNA index** | Bowtie2 index for rRNA contamination | Build or reuse existing |
| 8 | **PhiX index** | Bowtie2 index for PhiX spike-in | Reuse existing |

### Chromosome Naming Convention

Ensembl uses simple chromosome names (`1, 2, ..., 20, X, Y, MT`), but the pipeline expects UCSC-style `chr` prefixes (`chr1, chr2, ..., chr20, chrX, chrY, chrM`). The pipeline's QC step (`wdl/compute_mapped/mapped.wdl`) uses `grep "chrX"`, `grep "chrY"`, etc. to compute chromosome mapping statistics, so reference files **must** use chr-prefixed names.

A renaming script (`scripts/add_chr_prefix.sh`) is provided to convert Ensembl naming to UCSC convention after downloading. See Step 2 below.

---

## Current rn8 Reference Files

The rn8 (GRCr8, Ensembl 115) reference files are stored in GCS at:

```
gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/
```

| File | Size | Description |
|------|------|-------------|
| `Rattus_norvegicus.GRCr8.dna.toplevel.fa` | 2.70 GiB | Rat GRCr8 reference genome sequence (FASTA), with UCSC-style chr-prefixed chromosome names |
| `Rattus_norvegicus.GRCr8.115.gtf` | 586 MiB | Ensembl release 115 gene annotation for the rat GRCr8 genome, with UCSC-style chr-prefixed chromosome names |
| `refFlat_GRCr8_v115.txt` | 20.8 MiB | Gene intervals in refFlat format, used by Picard's CollectRnaSeqMetrics for RNA-seq QC |
| `rn8_v115_star_index.tar.gz` | 24.1 GiB | Pre-built STAR genome index for spliced read alignment |
| `rn8_rsem_reference.tar.gz` | 162 MiB | Pre-built RSEM reference index for transcript-level quantification (counts, TPM, FPKM) |

These files were built from Ensembl release 115, with chromosome names converted from Ensembl (`1, X, MT`) to UCSC (`chr1, chrX, chrM`) convention.

---

## Prerequisites

### Required Tools

#### 1. gtfToGenePred (refFlat generation)

Download from UCSC:
```bash
# Linux
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod +x gtfToGenePred
sudo mv gtfToGenePred /usr/local/bin/

# macOS (Intel)
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred
chmod +x gtfToGenePred
sudo mv gtfToGenePred /usr/local/bin/

# macOS (Apple Silicon) - use Rosetta or compile from source
# See: https://hgdownload.soe.ucsc.edu/admin/exe/
```

#### 2. Google Cloud SDK (gcloud storage)

Follow the official installation guide: https://cloud.google.com/sdk/docs/install

#### 3. Caper (WDL workflow runner)

Check the [caper repo](https://github.com/MoTrPAC/caper) for installation instructions.

**Verify all tools are installed:**
```bash
gtfToGenePred              # Should show usage
gcloud storage version             # Should show version
caper --version            # Should show version
```

---

## Step-by-Step Instructions

### Step 1: Download Ensembl Genome FASTA and GTF

Download the reference genome and gene annotations from the Ensembl FTP site. For rn8 (GRCr8, Ensembl release 115):

```bash
mkdir -p data/rat-ensembl-release-115 && cd data/rat-ensembl-release-115

# Download genome FASTA
wget https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz

# Download GTF annotation
wget https://ftp.ensembl.org/pub/release-115/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.115.gtf.gz

# Decompress
gunzip Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.GRCr8.115.gtf.gz

cd ..
```

**Verify the downloads:**
```bash
# Check chromosome names in FASTA (should be 1, 2, ..., X, Y, MT — bare Ensembl names)
grep '^>' data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.dna.toplevel.fa | head -5

# Check chromosome names in GTF (first column should be 1, 2, ..., X, Y, MT — bare Ensembl names)
grep -v '^#' data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.115.gtf | cut -f1 | sort -u | head
```

---

### Step 2: Add `chr` Prefix to Chromosome Names

Ensembl files use bare chromosome names (`1`, `X`, `MT`), but the pipeline requires UCSC-style names (`chr1`, `chrX`, `chrM`). Run the provided script to rename chromosomes in both the FASTA and GTF files:

```bash
bash scripts/add_chr_prefix.sh \
    data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.dna.toplevel.fa \
    data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.115.gtf
```

The script applies the following mapping (scaffolds/contigs are left unchanged):

| Ensembl | UCSC |
|---------|------|
| 1, 2, ..., 20 | chr1, chr2, ..., chr20 |
| X, Y | chrX, chrY |
| MT | chrM |

**Verify the renaming:**
```bash
# FASTA headers should now show chr1, chr2, etc.
grep '^>' data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.dna.toplevel.fa | head -5

# GTF first column should now show chr1, chr2, etc.
grep -v '^#' data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.115.gtf | cut -f1 | sort -uV

# Sanity check: no bare Ensembl chromosome names should remain
grep -v '^#' data/rat-ensembl-release-115/Rattus_norvegicus.GRCr8.115.gtf \
    | cut -f1 | grep -P '^([0-9]+|X|Y|MT)$' && echo "ERROR: bare names found" || echo "OK"
```

---

### Step 3: Generate refFlat File

The refFlat file is used by Picard's `CollectRnaSeqMetrics` for QC. Generate it from the GTF (which now has `chr` prefixes):

```bash
cd data/rat-ensembl-release-115/

# Convert GTF to genePred format (intermediate file)
gtfToGenePred -genePredExt Rattus_norvegicus.GRCr8.115.gtf GRCr8_genePred.txt

# Convert genePred to refFlat format
# refFlat has 11 columns: geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' \
    GRCr8_genePred.txt > refFlat_GRCr8_v115.txt

# Clean up intermediate file
rm GRCr8_genePred.txt

```

**Verify the output:**
```bash
# Check file format (should have 11 tab-separated columns)
head -3 refFlat_GRCr8_v115.txt | cut -f1-5

# Verify chromosome names (should be chr1, chr2, ..., chrX, chrY, chrM)
cut -f3 refFlat_GRCr8_v115.txt | sort -uV

# Count entries
wc -l refFlat_GRCr8_v115.txt
```

---

### Step 4: Upload Source Files to GCS

Upload the genome FASTA, GTF, and refFlat files to Google Cloud Storage (both in bash and fish):

**bash**:

```bash
GCS_PATH="gs://omicspipelines-public-resources/rnaseq/references/rat/rn8"

# Upload genome FASTA
gcloud storage cp data/Rattus_norvegicus.GRCr8.dna.toplevel.fa "$GCS_PATH/"

# Upload GTF annotation
gcloud storage cp data/Rattus_norvegicus.GRCr8.115.gtf "$GCS_PATH/"

# Upload refFlat file
gcloud storage cp data/refFlat_GRCr8_v115.txt "$GCS_PATH/"

# Verify uploads
gcloud storage ls "$GCS_PATH/"
```

**fish**:

```bash
# Set the variable for this session
set -g GCS_PATH "gs://omicspipelines-public-resources/rnaseq/references/rat/rn8"

# Upload genome FASTA
gcloud storage cp Rattus_norvegicus.GRCr8.dna.toplevel.fa "$GCS_PATH/"

# Upload GTF annotation
gcloud storage cp Rattus_norvegicus.GRCr8.115.gtf "$GCS_PATH/"

# Upload refFlat file
gcloud storage cp refFlat_GRCr8_v115.txt "$GCS_PATH/"

# Verify uploads
gcloud storage ls "$GCS_PATH/"
```

---

### Step 5: Build STAR Index

The STAR index is built using a WDL workflow. This runs on the cloud and requires significant resources (~120 GB RAM). The WDL reads the FASTA and GTF from GCS, so these must already be uploaded (Step 4).

**Input JSON file:** `examples/input_json/tasks/star_index_rn8_inputs.json`

```bash
caper submit wdl/star_ref/star_index.wdl \
    -i examples/input_json/tasks/star_index_rn8_inputs.json
```

After the workflow completes:
1. Locate the output tarball (`rn8_v115_star_index.tar.gz`)
2. Upload to GCS:
   ```bash
   gcloud storage cp /path/to/output/rn8_v115_star_index.tar.gz "$GCS_PATH/"
   ```

---

### Step 6: Build RSEM Reference

The RSEM reference is built using a WDL workflow. Like the STAR index, it reads from GCS.

**Input JSON file:** `examples/input_json/tasks/rsem_reference_rn8.json`

For this, you need to have [Caper](https://github.com/MoTrPAC/caper) installed and set up in a VM on the cloud. Then you clone this repo on the VM and run the following command:

```bash
caper submit motrpac-rna-seq-pipeline/wdl/rsem_index/rsem_reference.wdl \
    -i examples/input_json/tasks/rsem_reference_rn8.json
```

After the workflow completes:
1. Locate the output tarball (`rn8_rsem_reference.tar.gz`)
2. Upload to GCS:
   ```bash
   gcloud storage cp /path/to/output/rn8_rsem_reference.tar.gz "$GCS_PATH/"
   ```

---

### Step 7: Contamination Screening Indices

The pipeline requires Bowtie2 indices for globin, rRNA, and PhiX contamination screening.

For a **new genome version of an existing organism** (e.g., rn8 for rat), reuse the existing indices since the globin, rRNA, and PhiX sequences are conserved:

- Globin: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn_globin.tar.gz`
- rRNA: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn_rRNA.tar.gz`
- PhiX: `gs://omicspipelines-public-resources/rnaseq/references/rat/phix.tar.gz`

For a **new organism**, build these indices from scratch:

```bash
# Build bowtie2 index for globin sequences
caper submit wdl/bowtie2_index/bowtie2_index.wdl \
    -i examples/input_json/tasks/bowtie2_index_globin_inputs.json

# Build bowtie2 index for rRNA sequences
caper submit wdl/bowtie2_index/bowtie2_index.wdl \
    -i examples/input_json/tasks/bowtie2_index_rrna_inputs.json

# PhiX is organism-independent and can be reused
```

---

### Step 8: Update Pipeline Configuration

After all reference files are in GCS, add a version block to `scripts/make_json_rnaseq.py`:

1. Add the version to the `--version` argument choices
2. Add a reference mapping block in the `make_json_dict()` function

See the existing rn6, rn7, and rn8 blocks in the script for the pattern to follow.

---

## Verification Checklist

Before running the pipeline, verify all files are in place:

```bash
GCS_PATH="gs://omicspipelines-public-resources/rnaseq/references/rat/rn8"

# List all files
gcloud storage ls -l "$GCS_PATH/"

# Expected files:
# - Rattus_norvegicus.GRCr8.dna.toplevel.fa   (genome FASTA, chr-prefixed)
# - Rattus_norvegicus.GRCr8.115.gtf           (GTF annotation, chr-prefixed)
# - rn8_v115_star_index.tar.gz                (STAR index)
# - rn8_rsem_reference.tar.gz                 (RSEM reference)
# - refFlat_GRCr8_v115.txt                    (refFlat file, chr-prefixed)
```

**Verify chromosome naming:**
```bash
# Check FASTA headers (should be chr1, chr2, ..., chrX, chrY, chrM)
gcloud storage cat "$GCS_PATH"/Rattus_norvegicus.GRCr8.dna.toplevel.fa | head -1
# Should show: >chr1 dna:primary_assembly ...

# Check GTF chromosome column
gcloud storage cat "$GCS_PATH"/Rattus_norvegicus.GRCr8.115.gtf | grep -v '^#' | cut -f1 | sort -uV
# Should show: chr1, chr2, ..., chr20, chrM, chrX, chrY (plus scaffolds)
```

---

## Running the Pipeline

Once all reference files are ready, generate input JSON and run the pipeline:

```bash
python3 scripts/make_json_rnaseq.py \
    -g gs://your-bucket/fastq_raw \
    -o ./output \
    -r qc_report.csv \
    -a rat \
    -v rn8 \
    -n 1 \
    -p your-gcp-project \
    -d us-docker.pkg.dev/motrpac-portal/rnaseq
```

---

## Reference Genome Summary

All reference files stored in GCS use UCSC-style chr-prefixed chromosome names, regardless of the Ensembl source.

| Version | Assembly | Ensembl Release | Chromosome Names (in GCS) |
|---------|----------|-----------------|---------------------------|
| rn6 | Rnor_6.0 | 96 | chr1-chr20, chrX, chrY, chrM |
| rn7 | mRatBN7.2 | 108 | chr1-chr20, chrX, chrY, chrM |
| rn8 | GRCr8 | 115 | chr1-chr20, chrX, chrY, chrM |

---

## References

- Ensembl FTP: https://ftp.ensembl.org/
- UCSC Genome Browser tools: https://hgdownload.soe.ucsc.edu/admin/exe/
- Caper: https://github.com/MoTrPAC/caper
