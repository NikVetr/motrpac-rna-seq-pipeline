# Preparing Reference Files for the RNA-seq Pipeline

This document describes how to prepare reference genome files for running the MoTrPAC RNA-seq pipeline with a new organism or genome build.

The example below uses the rat rn8 (GRCr8) assembly, but the same process applies to any new reference genome.

---

## Overview

The RNA-seq pipeline requires **7 reference files**, all using **UCSC-style chromosome naming** (chr1, chrX, etc.):

| # | File | Description | How to Create |
|---|------|-------------|---------------|
| 1 | **GTF file** | Gene annotation in GTF format (UCSC naming) | Convert from GFF3, then rename chromosomes |
| 2 | **STAR index** | Genome index for STAR aligner | Build with WDL workflow |
| 3 | **RSEM reference** | Reference for RSEM quantification | Build with WDL workflow |
| 4 | **refFlat file** | Gene intervals for Picard QC | Generate from GTF |
| 5 | **Globin index** | Bowtie2 index for globin contamination | Build or reuse existing |
| 6 | **rRNA index** | Bowtie2 index for rRNA contamination | Build or reuse existing |
| 7 | **PhiX index** | Bowtie2 index for PhiX spike-in | Reuse existing |

> ⚠️ **Important:** NCBI genomes use RefSeq accession-based chromosome names (e.g., `NC_086019.1`). These must be converted to UCSC-style names (e.g., `chr1`) before use in this pipeline. See Step 2 below.

---

## Prerequisites

### Required Tools

Install the following tools before starting:

#### 1. gffread (GFF3 to GTF conversion)

Download the pre-compiled binary:
```bash
# Download latest release (check https://github.com/gpertea/gffread/releases for current version)
curl -L https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz -o gffread.tar.gz
tar -xzf gffread.tar.gz
sudo mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin/

# For macOS, compile from source:
git clone https://github.com/gpertea/gffread.git
cd gffread
make release
sudo mv gffread /usr/local/bin/
```

Or install via Homebrew (preferred on macOS):
```bash
brew install gffread
```

#### 2. gtfToGenePred (refFlat generation)

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

#### 3. Google Cloud SDK (gsutil)

Follow the official installation guide: https://cloud.google.com/sdk/docs/install

```bash
# Verify installation
gsutil version
```

#### 4. Caper (WDL workflow runner): on a VM

Check the [caper repo](https://github.com/MoTrPAC/caper) to find out more about the installations


**Verify all tools are installed:**
```bash
gffread --version          # Should show version
gtfToGenePred              # Should show usage
gsutil version             # Should show version
caper --version            # Should show version
```

### Required Input Files

You need:
1. **Genomic FASTA file** - The reference genome sequence (`.fna` or `.fa`)
2. **Annotation file** - Gene annotations in GFF3 or GTF format

For NCBI genomes, download from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets).

---

## Step-by-Step Instructions

Download the NCBI genome. For example, to download the rat RN8

```
# create a "data" folder in this repo, `cd data` and download the genome like this:

datasets download genome accession GCF_036323735.1 --include gff3,rna,cds,protein,genome,seq-report
```

### Step 1: Convert GFF3 to GTF (if needed)

The pipeline requires **GTF format**. If your annotation is in GFF3 format, convert it:

```bash
# Convert GFF3 to GTF using gffread
gffread data/GCF_036323735.1/genomic.gff \
    -T \
    -o data/GCF_036323735.1/Rattus_norvegicus.GRCr8.gtf
```

**Verify the output:**
```bash
# Check file was created and has content
wc -l data/GCF_036323735.1/Rattus_norvegicus.GRCr8.gtf

# Preview the GTF format (should have gene_id and transcript_id attributes)
head -5 data/GCF_036323735.1/Rattus_norvegicus.GRCr8.gtf
```

---

### Step 2: Rename Chromosomes to UCSC-style Naming (CRITICAL)

> ⚠️ **This step is required for NCBI genomes.** NCBI uses RefSeq accession-based chromosome names (e.g., `NC_086019.1`), but the RNA-seq pipeline expects UCSC-style names (e.g., `chr1`, `chrX`, `chrM`). Skipping this step will cause pipeline failures.

The `sequence_report.jsonl` file from NCBI Datasets contains the mapping between RefSeq accessions and UCSC-style names. We provide a script to automate this conversion.

#### Run the chromosome renaming script:

```bash
# Create output directory for UCSC-named files
mkdir -p data/GCF_036323735.1_ucsc

# Run the renaming script
python3 scripts/rename_chromosomes.py \
    data/GCF_036323735.1/sequence_report.jsonl \
    data/GCF_036323735.1 \
    data/GCF_036323735.1_ucsc
```

The script will:
1. Parse the chromosome mapping from `sequence_report.jsonl`
2. Rename chromosomes in the genomic FASTA file
3. Rename chromosomes in the GTF annotation file
4. Perform QC validation checks
5. Generate a summary report and genes-per-chromosome plot

#### Expected output:

```
data/GCF_036323735.1_ucsc/
├── GCF_036323735.1_GRCr8_genomic_ucsc.fna   # Renamed FASTA
├── Rattus_norvegicus.GRCr8_ucsc.gtf         # Renamed GTF
├── chromosome_mapping.tsv                    # Mapping file for reference
├── qc_report.txt                             # QC validation report
└── genes_per_chromosome.png                  # QC plot
```

#### QC Validation

The script performs automatic validation:

- ✓ **Sequence count**: Verifies the same number of sequences in input/output FASTA
- ✓ **File size**: Confirms file sizes match within 1% (small differences from header changes)
- ✓ **Main chromosomes**: Checks that chr1, chr2, chrX, chrY exist in output
- ✓ **Gene distribution**: Generates a plot showing genes per chromosome

**Review the QC report before proceeding:**

```bash
cat data/GCF_036323735.1_ucsc/qc_report.txt
```

Sample output:
```
============================================================
CHROMOSOME RENAMING QC REPORT
============================================================

FASTA File:
  Total sequences: 76
  Renamed sequences: 76
  Input file size: 2,841,234,567 bytes
  Output file size: 2,841,234,123 bytes
  ✓ File size difference: 0.001% (OK)
  ✓ All sequences were mapped

GTF File:
  Total lines: 1,234,567
  Renamed lines: 1,234,567
  ✓ All chromosomes were mapped

Genes per Chromosome:
  chr1                       3,245 ████████████
  chr2                       2,891 ██████████
  ...
  chrX                       1,234 █████
  chrY                         456 ██

Total genes: 47,357
============================================================
```

**Important:** Use the files from `data/GCF_036323735.1_ucsc/` for all subsequent steps.

---

### Step 3: Generate refFlat File

The refFlat file is used by Picard's `CollectRnaSeqMetrics` for QC. Generate it from the **UCSC-renamed GTF**:

```bash
# Navigate to the UCSC-renamed data directory
cd data/GCF_036323735.1_ucsc

# Convert GTF to genePred format (intermediate file)
gtfToGenePred -genePredExt Rattus_norvegicus.GRCr8_ucsc.gtf rn8_genePred.txt

# Convert genePred to refFlat format
# refFlat has 11 columns: geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' \
    rn8_genePred.txt > refFlat_rn8_GRCr8.txt

# Clean up intermediate file
rm rn8_genePred.txt

# Return to repo root
cd ../..
```

**Verify the output:**

```bash
# Check file format (should have 11 tab-separated columns with chr-prefixed names)
head -3 data/GCF_036323735.1_ucsc/refFlat_rn8_GRCr8.txt | cut -f1-5

# Verify chromosomes have UCSC naming
cut -f3 data/GCF_036323735.1_ucsc/refFlat_rn8_GRCr8.txt | sort -u | head -10

# Count entries
wc -l data/GCF_036323735.1_ucsc/refFlat_rn8_GRCr8.txt
```

---

### Step 3: Upload Source Files to GCS

Upload the **UCSC-renamed** genomic FASTA, GTF, and refFlat files to Google Cloud Storage:

```bash
# Set the destination path
GCS_PATH="gs://omicspipelines-public-resources/rnaseq/references/rat/rn8"

# Upload UCSC-renamed genomic FASTA
gsutil cp data/GCF_036323735.1_ucsc/GCF_036323735.1_GRCr8_genomic_ucsc.fna ${GCS_PATH}/

# Upload UCSC-renamed GTF annotation
gsutil cp data/GCF_036323735.1_ucsc/Rattus_norvegicus.GRCr8_ucsc.gtf ${GCS_PATH}/

# Upload refFlat file (generated from UCSC-renamed GTF)
gsutil cp data/GCF_036323735.1_ucsc/refFlat_rn8_GRCr8.txt ${GCS_PATH}/

# Upload chromosome mapping for reference
gsutil cp data/GCF_036323735.1_ucsc/chromosome_mapping.tsv ${GCS_PATH}/
```

**Verify uploads:**
```bash
gsutil ls ${GCS_PATH}/
```

---

### Step 5: Build STAR Index

The STAR index is built using a WDL workflow. This runs on the cloud and requires significant resources (~120 GB RAM).

**Input JSON file:** `examples/input_json/tasks/star_index_rn8_inputs.json`

```bash
# Run the STAR index workflow with Caper
caper submit wdl/star_ref/star_index.wdl \
    -i examples/input_json/tasks/star_index_rn8_inputs.json
```

After the workflow completes:
1. Locate the output tarball (`rn8_GRCr8_star_index.tar.gz`)
2. Upload to GCS:
   ```bash
   gsutil cp /path/to/output/rn8_GRCr8_star_index.tar.gz ${GCS_PATH}/
   ```

---

### Step 6: Build RSEM Reference

The RSEM reference is built using a WDL workflow.

**Input JSON file:** `examples/input_json/tasks/rsem_reference_rn8.json`

```bash
# Run the RSEM reference workflow with Caper
caper submit wdl/rsem_index/rsem_reference.wdl \
    -i examples/input_json/tasks/rsem_reference_rn8.json
```

After the workflow completes:
1. Locate the output tarball (`rn8_rsem_reference.tar.gz`)
2. Upload to GCS:
   ```bash
   gsutil cp /path/to/output/rn8_rsem_reference.tar.gz ${GCS_PATH}/
   ```

---

### Step 7: Contamination Screening Indices

The pipeline requires bowtie2 indices for globin, rRNA, and PhiX contamination screening. These are **required** inputs.

For a **new genome version of an existing organism** (e.g., rn8 for rat), you can reuse the existing indices since the globin, rRNA, and PhiX sequences are conserved:

- Globin: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn_globin.tar.gz`
- rRNA: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn_rRNA.tar.gz`
- PhiX: `gs://omicspipelines-public-resources/rnaseq/references/rat/phix.tar.gz`

For a **new organism**, you must build these indices from scratch:

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

After all reference files are in GCS, update `scripts/make_json_rnaseq.py` to add support for the new genome version:

1. Add the version to the `--version` argument choices
2. Add a reference mapping block in the `make_json_dict()` function

See the existing rn6, rn7, and rn8 blocks in the script for examples.

---

## Verification Checklist

Before running the pipeline, verify all files are in place:

```bash
GCS_PATH="gs://omicspipelines-public-resources/rnaseq/references/rat/rn8"

# List all files
gsutil ls -l ${GCS_PATH}/

# Expected files (note: using UCSC-style chromosome naming):
# - GCF_036323735.1_GRCr8_genomic_ucsc.fna   (FASTA with chr1, chrX, etc.)
# - Rattus_norvegicus.GRCr8_ucsc.gtf         (GTF with chr1, chrX, etc.)
# - rn8_GRCr8_star_index.tar.gz              (STAR index)
# - rn8_rsem_reference.tar.gz                (RSEM reference)
# - refFlat_rn8_GRCr8.txt                    (refFlat file)
# - chromosome_mapping.tsv                    (RefSeq->UCSC name mapping)
```

**Verify chromosome naming is correct:**
```bash
# Check FASTA headers
gsutil cat ${GCS_PATH}/GCF_036323735.1_GRCr8_genomic_ucsc.fna | head -1
# Should show: >chr1 ...

# Check GTF chromosome column
gsutil cat ${GCS_PATH}/Rattus_norvegicus.GRCr8_ucsc.gtf | head -5
# First column should be chr1, chr2, etc.
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

## rn8 Assembly Information

This directory contains the rat rn8 (GRCr8) reference genome:

- **Assembly Name:** GRCr8 (GCF_036323735.1)
- **Organism:** Rattus norvegicus (Norway rat)
- **Strain:** BN/NHsdMcwi (Brown Norway)
- **Release Date:** January 31, 2024
- **Annotation:** NCBI RefSeq (GCF_036323735.1-RS_2024_02)
- **Gene Count:** 47,357 total genes (23,154 protein-coding)

### Directory Contents

If you download the data from the NCBI as explained above, this is the typical content:

```
data/
├── README.md                              # This file
├── assembly_data_report.jsonl             # Assembly metadata
├── dataset_catalog.json                   # File manifest
└── GCF_036323735.1/
    ├── GCF_036323735.1_GRCr8_genomic.fna  # Genomic sequence (FASTA)
    ├── genomic.gff                        # Gene annotation (GFF3, original)
    ├── Rattus_norvegicus.GRCr8.gtf        # Gene annotation (GTF, converted)
    ├── refFlat_rn8_GRCr8.txt              # refFlat file (generated)
    ├── cds_from_genomic.fna               # CDS sequences
    ├── rna.fna                            # RNA sequences
    ├── protein.faa                        # Protein sequences
    └── sequence_report.jsonl              # Sequence metadata
```

---


## References

- NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets
- gffread: https://github.com/gpertea/gffread
- UCSC Genome Browser tools: https://hgdownload.soe.ucsc.edu/admin/exe/
- Caper: https://github.com/MoTrPAC/caper
