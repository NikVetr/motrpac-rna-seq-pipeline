MoTrPAC RNA-SEQ Pipeline
=================================================

[![DOI](https://zenodo.org/badge/144631622.svg)](https://zenodo.org/badge/latestdoi/144631622)

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Prerequisites](#prerequisites)
- [GCP Set-up](#gcp-set-up)
- [Software / Dockerfiles](#software--dockerfiles)
- [Configuration Files](#configuration-files)
- [Run the Pipeline](#run-the-pipeline)
- [Pipeline Outputs](#pipeline-outputs)
- [Monitoring and Job Management](#monitoring-and-job-management)
- [WDL Workflow Structure](#wdl-workflow-structure)
- [Local Development and Testing](#local-development-and-testing)
- [Troubleshooting](#troubleshooting)
- [Citations and References](#citations-and-references)
- [Contributing and Support](#contributing-and-support)
- [Version Information](#version-information)

## Overview

This repo contains the rna-seq data processing pipeline implemented in Workflow Description Language (WDL) based on harmonized [RNA-SEQ MOP](https://docs.google.com/document/d/e/2PACX-1vRFurZraZfxfMd5BWfIQEnETlalDNjQPyMjS7TCTgc3MMlMtB_-tmJfEK7lmRV7GD30I7R9-ISX3kuM/pub). This pipeline uses [caper](https://github.com/MoTrPAC/caper), a wrapper python package for the workflow management system [Cromwell](https://cromwell.readthedocs.io/en/stable/). All the data was processed on the Google Cloud Platform (GCP).

### Supported Organisms and Genome Builds

The pipeline supports the following organisms and genome versions:
- **Rat**: rn6 (Rnor_6.0, Ensembl 96), rn7 (mRatBN7.2, Ensembl 108), rn8 (GRCr8, Ensembl 115)
- **Human**: GENCODE v39 or v47 (GRCh38)

### Pipeline Tools

The pipeline uses:
- [STAR aligner](https://github.com/alexdobin/STAR) for read alignment
- [RSEM](https://github.com/deweylab/RSEM) for transcript quantification (TPM, FPKM, counts)
- [featureCounts](https://subread.sourceforge.net/) for gene-level read quantification
- [UMI-tools](https://umi-tools.readthedocs.io/) for directional UMI grouping
- Quality control tools including FastQC, MultiQC, Picard, and Bowtie2 (for contamination assessment)

### Pipeline Outputs

The pipeline generates:
- Gene expression quantification (counts, TPM, FPKM) from RSEM
- Gene counts from featureCounts
- Comprehensive QC metrics for outlier detection and covariate adjustment
- Optional legacy-compatible MultiQC reports for pre- and post-alignment QC

## Quick Start

For experienced users, here's the essential workflow:

```bash
# 1. Clone the repository
git clone https://github.com/MoTrPAC/motrpac-rna-seq-pipeline

# 2. Install Python dependencies
pip3 install -r scripts/requirements.txt

# 3. Generate input JSON configuration
python3 scripts/make_json_rnaseq.py \
  -g gs://your-bucket/fastq_raw \
  -o ./input_json \
  -r batch1_qc_metrics \
  -a human \
  -v gencode_v47 \
  -n 1 \
  -p your-gcp-project \
  -d us-docker.pkg.dev/motrpac-portal/rnaseq

# 4. Submit the pipeline
caper submit wdl/rnaseq_pipeline_scatter.wdl -i input_json/set1_rnaseq.json

# 5. Monitor pipeline status
caper list
```

## Prerequisites

### Required Accounts
- Google Cloud Platform (GCP) account with billing enabled
- GCP service account with appropriate permissions
- GCP Storage bucket for pipeline inputs and outputs

### Required Software (Local Machine)
- [Google Cloud SDK](https://cloud.google.com/sdk/docs/install)
- Python >= 3.6.9
- Git

### Python Dependencies
Install required Python packages:
```bash
pip3 install -r scripts/requirements.txt
```

The main dependencies include:
- `gcsfs` - for accessing Google Cloud Storage
- `numpy` - for data processing

### GCP Permissions and APIs
Ensure the following APIs are enabled in your GCP project:
- Compute Engine API
- Cloud Storage API
- Google Cloud Batch API (for workflow execution)

## GCP Set-up

The WDL/Cromwell framework is optimized to run pipelines in high-performance computing environments. The MoTrPAC Bioinformatics Center runs pipelines on Google Cloud Platform (GCP). We used a number of fantastic tools developed by our colleagues from the [ENCODE project](https://github.com/ENCODE-DCC) to run pipelines on GCP (and other HPC platforms).

A brief summary of the steps to set-up a VM to run the Motrpac pipelines on GCP (**for details, please, check the [caper repo](https://github.com/MoTrPAC/caper/blob/master/scripts/gcp_caper_server/README.md)**):

- Create a GCP account.
- Enable cloud APIs.
- Install the [Google Cloud SDK](https://cloud.google.com/sdk/docs/install) (Software Development Kit) on your local machine.
- Create a service account and download the key file to your local computer (e.g. `service-account-191919.json`)
- Create a bucket for pipeline inputs and outputs (e.g. gs://pipelines/). Note: a GCP bucket is similar to a folder on your computer or a storage unit, but it is stored on Google's servers in the cloud instead of on your local computer.
- Set up a VM on GCP: create a Virtual Machine (VM) instance from where the pipelines will be run. We recommend the script available in the [caper repo](https://github.com/MoTrPAC/caper/). For that, clone the repo on your local machine and run the following command:

 ```bash
 $ bash create_instance.sh [INSTANCE_NAME] [PROJECT_ID] [GCP_SERVICE_ACCOUNT_KEY_JSON_FILE] [GCP_OUT_DIR]

 # Example for the pipeline:
./create_instance.sh pipeline-instance your-gcp-project-name service-account-191919.json gs://pipelines/results/
```

- Finally, clone the repo on your VM instance

 ```bash
 git clone https://github.com/MoTrPAC/motrpac-rna-seq-pipeline
 ```

## Software / Dockerfiles

Several tools are required to run the rna-seq pipeline. All of them are pre-installed in docker containers, which are publicly available in the [Artifact Registry](https://cloud.google.com/artifact-registry).

### Available Docker Images

The pipeline uses the following containerized tools (all available at `us-docker.pkg.dev/motrpac-portal/rnaseq`):

- `fastqc:latest` - FastQC for quality control
- `umi_attach:latest` - UMI attachment utility
- `cutadapt:latest` - Adapter trimming
- `multiqc:latest` - Aggregate QC reporting
- `star:latest` - STAR aligner
- `feature_counts:latest` - featureCounts from Subread
- `rsem:latest` - RSEM quantification
- `bowtie:latest` - Bowtie2 aligner (for contamination screening)
- `picard:latest` - Picard tools (MarkDuplicates, CollectRnaSeqMetrics)
- `umi_dup:latest` - UMI-based duplication assessment
- `samtools:latest` - SAMtools utilities
- `collect_qc:latest` - Custom QC metrics collection
- `merge_results:latest` - Result merging across samples

### Building and Updating Containers

To find out more about the specific versions of tools used to run the pipeline, check the `dockerfiles/*.Dockerfile`.

To build and push updated containers:
```bash
# Build all dockerfiles
bash scripts/build_dockerfiles.sh

# Push to Artifact Registry (requires appropriate permissions)
bash scripts/push_dockerfiles.sh
```

## Configuration Files

An input configuration file (in JSON format) is required to process the data through the pipeline. This configuration file contains several key-value pairs that specify the inputs and outputs of the workflow, the location of the input files, default pipeline parameters, docker containers, the execution environment, and other parameters needed for execution.

### Generating Configuration Files

The optimal way to generate the configuration files is to run the `make_json_rnaseq.py` script.

**Usage:**
```bash
python3 scripts/make_json_rnaseq.py \
  -g GCP_PATH \               # GCS path to directory containing FASTQ files
  -o OUTPUT_PATH \            # Local path where JSON files will be written
  -r OUTPUT_REPORT_NAME \     # Name for the output QC metrics report
  -a {rat,human} \            # Organism
  -v {rn6,rn7,rn8,gencode_v39,gencode_v47} \  # Genome/annotation version
  -n NUM_CHUNKS \             # Number of batches to split samples into
  -p PROJECT \                # GCP project name
  -d DOCKER_REPO \            # Docker repository prefix (optional)
  -u                          # Include undetermined reads (optional)
```

**Complete Example:**
```bash
python3 scripts/make_json_rnaseq.py \
  -g gs://motrpac-bucket/rna-seq/human/batch7_20220316/fastq_raw \
  -o ./input_json \
  -r batch7_qc_metrics.csv \
  -a human \
  -v gencode_v47 \
  -n 1 \
  -p motrpac-portal \
  -d us-docker.pkg.dev/motrpac-portal/rnaseq \
  -i
```

This will create JSON configuration file(s) (e.g., `set1_rnaseq.json`, `set2_rnaseq.json`, etc.) in the specified output directory.

### Modernization controls

GENCODE v47 automatically selects the immutable release profile in
`config/release-profiles/human-gencode-v47.json`; v39 retains the historical
references and images. Matched I1 reads and directional UMI molecule-expression
matrices are enabled by default, while the conventional all-read RSEM and
featureCounts outputs are retained. Pass `--legacy-all-read-expression-only`
to omit only the molecule-level matrices.

Pre-trim FastQC, post-trim FastQC, contamination QC, alignment QC, and UMI QC
remain enabled by default and can be disabled independently with the documented
`--skip-*` options. `--combine-contamination-qc` runs the three screens on one
worker; adding `--contamination-qc-pairs N` deterministically samples the same
post-trim read pairs for all three screens (`N=0` retains full depth).
The QC table is assembled directly from native tool reports. Pass
`--run-multiqc` to additionally publish the legacy pre- and post-alignment
MultiQC archives; this requires both FastQC groups and alignment QC.

For GCP tests, select a complete profile with `--runtime-profile` and select
STAR scratch explicitly with `--star-disk-type HDD|SSD`. These options augment
the existing JSON format and Caper submission command; see
[`scripts/scripts_readme.md`](scripts/scripts_readme.md) for the full generator
interface and [`docs/gcp-cli-canary-runbook.md`](docs/gcp-cli-canary-runbook.md)
for the bounded benchmark procedure.

### Organism Reference Files

The `make_json_rnaseq.py` script automatically selects the appropriate reference files based on the organism and version:

**Rat (rn6):**
- STAR index: `gs://omicspipelines-public-resources/rnaseq/references/rat/Rnor6_v96_star_index.tar.gz`
- GTF: `gs://omicspipelines-public-resources/rnaseq/references/rat/Rattus_norvegicus.Rnor_6.0.96.gtf`
- RSEM reference: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn6_rsem_reference.tar.gz`

**Rat (rn7):**
- STAR index: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/rn7_v108_star_index.tar.gz`
- GTF: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/Rattus_norvegicus.mRatBN7.2.108.gtf`
- RSEM reference: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn7/rn7_rsem_reference.tar.gz`

**Rat (rn8):**
- STAR index: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/rn8_v115_star_index.tar.gz`
- GTF: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/Rattus_norvegicus.GRCr8.115.gtf`
- RSEM reference: `gs://omicspipelines-public-resources/rnaseq/references/rat/rn8/rn8_rsem_reference.tar.gz`

**Human (gencode_v39):**
- STAR index: `gs://omicspipelines-public-resources/rnaseq/references/human/hg38_v39_star_index.tar.gz`
- GTF: `gs://omicspipelines-public-resources/rnaseq/references/human/GRCh38.v39.primary_assembly.annotation.gtf`
- RSEM reference: `gs://omicspipelines-public-resources/rnaseq/references/human/hg38_rsem_reference.tar.gz`

**Human (gencode_v47):**
- Immutable references and matching tool images are defined together in
  `config/release-profiles/human-gencode-v47.json`.

For more details, see the [scripts documentation](scripts/scripts_readme.md).

## Run the Pipeline

Connect to the VM and submit the job using the below command:

```bash
caper submit wdl/rnaseq_pipeline_scatter.wdl -i input_json/set1_rnaseq.json
```

Check the status of workflows and make sure they have succeeded by typing `caper list` on the VM instance that's running the job and look for `Succeeded`.

## Pipeline Outputs

The pipeline generates the following main output files:

### Quantification Files

1. **RSEM Gene Expression Quantification**
   - `*_rsem_genes_count.txt` - Raw gene-level counts
   - `*_rsem_genes_tpm.txt` - Transcripts Per Million (TPM) normalized expression
   - `*_rsem_genes_fpkm.txt` - Fragments Per Kilobase Million (FPKM) normalized expression

2. **featureCounts Gene Quantification**
   - `*_feature_counts.txt` - Gene-level raw counts from featureCounts

### Quality Control Files

3. **QC Metrics Report**
   - `*_qc_report.csv` - Comprehensive QC metrics per sample including:
     - Read alignment statistics
     - rRNA, globin, and PhiX contamination rates
     - PCR duplication rates
     - Strand specificity
     - 5' to 3' coverage bias
     - Percentage of reads mapping to coding/intronic/intergenic regions
     - Chromosome mapping percentages

### Additional Outputs (Per Sample)

The pipeline also generates intermediate outputs for each sample (stored in Cromwell execution directories):
- FastQC reports (pre- and post-trimming)
- STAR alignment BAM files
- Trimmed FASTQ files
- Picard metrics files

When `--run-multiqc` is selected, the pre- and post-alignment consolidated
report archives are also published as top-level outputs.

### Retrieving Outputs

Final merged outputs are written to the GCS bucket specified during pipeline submission. Individual sample outputs are organized in the Cromwell execution directory structure.

## Monitoring and Job Management

### Checking Pipeline Status

```bash
# List all workflows
caper list

# Check detailed status of a specific workflow
caper metadata [WORKFLOW_ID]
```

### Monitoring Running Jobs

```bash
# View workflows currently running
caper list | grep Running

# Check logs for a specific workflow
caper debug [WORKFLOW_ID]
```

### Managing Workflows

```bash
# Abort a running workflow
caper abort [WORKFLOW_ID]

# Check troubleshooting information
caper troubleshoot [WORKFLOW_ID]
```

### Retrieving Results

Successful pipeline runs will write outputs to your specified GCS bucket. Intermediate files and execution logs are stored in:
- `cromwell-executions/` - Contains all task execution outputs and logs
- `cromwell-workflow-logs/` - Contains workflow-level logs

To copy results from GCS to your local machine:
```bash
gsutil -m cp -r gs://your-bucket/results/workflow_id/* ./local_results/
```

## WDL Workflow Structure

The pipeline is organized as a modular WDL workflow with the following structure:

### Main Workflow
- `wdl/rnaseq_pipeline_scatter.wdl` - Main workflow that orchestrates all tasks using a scatter-gather pattern to process multiple samples in parallel

### Task Modules

The pipeline consists of the following task modules (in `wdl/` directory):

**Pre-alignment QC and Processing:**
- `fastqc/` - Quality control with FastQC (pre- and post-trimming)
- `attach_umi/` - Attach UMI indices to read names
- `cutadapt/` - Adapter trimming
- `multiqc/` - Aggregate QC reporting

**Alignment and Quantification:**
- `star_align/` - Alignment with STAR
- `rsem_exp/` - RSEM quantification
- `feature_counts/` - featureCounts quantification

**Contamination Screening:**
- `bowtie2_align/` - Bowtie2 alignment to globin, rRNA, and PhiX references

**Post-alignment QC:**
- `mark_duplicates/` - PCR duplicate marking with Picard
- `collect_rnaseq_metrics/` - RNA-seq QC metrics with Picard
- `umi_dup/` - Directional UMI grouping, duplication QC, and molecule BAM preparation
- `compute_mapped/` - Chromosome mapping statistics
- `collect_qc_metrics/` - Consolidated QC metrics collection

**Results Aggregation:**
- `merge_results/` - Merge quantification and QC data across all samples

### Reference Building Workflows

The repository also includes workflows for building reference files:
- `wdl/star_ref/` - Build STAR genome indices
- `wdl/rsem_index/` - Build RSEM reference indices
- `wdl/bowtie2_index/` - Build Bowtie2 indices

For detailed instructions on how to prepare all required reference files for a new organism or genome build (index building, refFlat generation, etc.), see [README-DATA-REF.md](README-DATA-REF.md).

## Local Development and Testing

### Setup for Development

Use the provided setup scripts to configure your development environment:

```bash
# Set up VM for pipeline execution
bash scripts/setup/setup_vm.sh

# Set up local development environment
bash scripts/setup/setup_develop.sh
```

### Validating JSON Files

Before submitting workflows, validate your JSON configuration files:

```bash
python3 scripts/validate_jsons.py input_json/set1_rnaseq.json
```

### Testing with Prototype Examples

The `prototype/` directory contains example configuration files and submission scripts:

```bash
# Example submission script for generic use
bash prototype/submit_rnaseq_generic.sh

# Example JSON configurations in prototype/input_json/
```

The `examples/` directory contains additional JSON examples for individual tasks and different organism configurations.

### Building Docker Images Locally

```bash
# Build all docker images
bash scripts/build_dockerfiles.sh

# Push to your container registry (configure registry URL first)
bash scripts/push_dockerfiles.sh
```

## Troubleshooting

### Common Issues

**1. Pipeline Fails During Submission**
- Verify JSON configuration is valid using `scripts/validate_jsons.py`
- Ensure all required input files exist in the specified GCS paths
- Check that service account has permissions to access GCS buckets

**2. Tasks Fail with "Out of Memory" Errors**
- Increase `*_ramGB` parameters in your JSON configuration
- Default memory allocations are in `scripts/make_json_rnaseq.py`

**3. Tasks Fail with "Out of Disk Space" Errors**
- Increase `*_disk` parameters in your JSON configuration
- Ensure your GCS bucket has sufficient quota

**4. Cannot Find Output Files**
- Check workflow succeeded: `caper list`
- Outputs are in the GCS bucket specified in your configuration
- Check Cromwell execution logs in `cromwell-executions/`

**5. Docker Image Pull Failures**
- Verify you have access to the Artifact Registry
- Check that docker image names/tags are correct in JSON
- Ensure Compute Engine service account has Artifact Registry Reader role

### Accessing Logs

**Workflow-level logs:**
```bash
# View workflow metadata
caper metadata [WORKFLOW_ID]

# Check troubleshooting info
caper troubleshoot [WORKFLOW_ID]
```

**Task-level logs:**
Navigate to the Cromwell execution directory:
```bash
cd cromwell-executions/rnaseq_pipeline/[WORKFLOW_ID]/
# Find specific task directories and check stderr/stdout logs
```

**GCP Console:**
- Navigate to Life Sciences API in GCP Console
- View operation logs and details for each task execution

### Getting Help

If issues persist:
1. Check the Cromwell documentation: https://cromwell.readthedocs.io/
2. Review the Caper documentation: https://github.com/MoTrPAC/caper/
3. Open an issue on the GitHub repository with:
   - Workflow ID
   - Error messages from logs
   - JSON configuration (with sensitive data removed)

## Citations and References

### Pipeline Documentation
- [RNA-SEQ MOP](https://docs.google.com/document/d/e/2PACX-1vRFurZraZfxfMd5BWfIQEnETlalDNjQPyMjS7TCTgc3MMlMtB_-tmJfEK7lmRV7GD30I7R9-ISX3kuM/pub) - MoTrPAC RNA-seq Method of Procedure

### Workflow Management
- [Cromwell](https://cromwell.readthedocs.io/en/stable/) - Workflow management system
- [Caper](https://github.com/MoTrPAC/caper/) - Cromwell wrapper for easy workflow execution
- [WDL](https://openwdl.org/) - Workflow Description Language specification

### Analysis Tools
- [STAR](https://github.com/alexdobin/STAR) - Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013.
- [RSEM](https://github.com/deweylab/RSEM) - Li B and Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011.
- [featureCounts](https://subread.sourceforge.net/) - Liao Y, et al. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control for high throughput sequence data
- [MultiQC](https://multiqc.info/) - Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016.
- [Cutadapt](https://cutadapt.readthedocs.io/) - Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011.
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) - Langmead B and Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012.
- [Picard Tools](https://broadinstitute.github.io/picard/) - Broad Institute toolkit for SAM/BAM file manipulation

### Infrastructure
- [ENCODE-DCC](https://github.com/ENCODE-DCC) - Tools and pipelines from the ENCODE Project Consortium

### Reference Genomes
- **Rat rn6**: Ensembl Rnor_6.0 release 96
- **Rat rn7**: Ensembl mRatBN7.2 release 108
- **Rat rn8**: Ensembl GRCr8 release 115
- **Human**: GENCODE v39 or v47 (GRCh38)

## Contributing and Support

### Reporting Issues

If you encounter bugs or have feature requests, please open an issue on the [GitHub repository](https://github.com/MoTrPAC/motrpac-rna-seq-pipeline/issues).

When reporting issues, please include:
- Description of the problem
- Steps to reproduce
- Expected vs. actual behavior
- Workflow ID (if applicable)
- Relevant error messages or logs
- JSON configuration (remove sensitive information)

### Contact

For questions or support related to the MoTrPAC RNA-seq pipeline:
- Open an issue on GitHub: https://github.com/MoTrPAC/motrpac-rna-seq-pipeline/issues
- Contact the MoTrPAC Bioinformatics Center

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request with a clear description of the changes

### Related Repositories

- [MoTrPAC Data Hub](https://motrpac-data.org/) - Access to MoTrPAC datasets
- [MoTrPAC GitHub Organization](https://github.com/MoTrPAC) - Other MoTrPAC analysis pipelines and tools

## Version Information

### Current Version
This pipeline is actively maintained and updated. Check the [releases page](https://github.com/MoTrPAC/motrpac-rna-seq-pipeline/releases) for version history and changelogs.

### Citing This Pipeline

If you use this pipeline in your research, please cite:

[![DOI](https://zenodo.org/badge/144631622.svg)](https://zenodo.org/badge/latestdoi/144631622)

### Compatibility Notes

- **WDL Version**: 1.0
- **Cromwell Version**: Compatible with Cromwell 50+
- **Python Version**: Requires Python >= 3.6.9
- **GCP**: Designed for Google Cloud Platform (adaptable to other backends with Cromwell configuration)

### Change History

Major updates and changes are documented in the repository's commit history. For significant changes:
- Reference genome updates (rn6 → rn7, GENCODE versions)
- Tool version updates (see dockerfiles for current versions)
- Workflow optimizations and bug fixes

Check the [commit history](https://github.com/MoTrPAC/motrpac-rna-seq-pipeline/commits/) for detailed changes.
