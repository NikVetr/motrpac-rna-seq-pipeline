# Scripts

### `make_json_rnaseq.py`

Generates the input configuration file required to run the rna-seq pipeline.

- Requires Python `>3.6.9`
- Install required packages by running `pip3 install -r scripts/requirements.txt`

```
usage: make_json_rnaseq.py [-h] -g GCP_PATH -o OUTPUT_PATH
                           -r OUTPUT_REPORT_NAME [-u] -a {rat,human}
                           -v {rn6,rn7,rn8,gencode_v39,gencode_v47}
                           -n NUM_CHUNKS [-d DOCKER_REPO]
                           [--release-manifest RELEASE_MANIFEST]
                           [--runtime-profile RUNTIME_PROFILE]
                           [--star-disk-type {HDD,SSD}] [-i]
                           [--legacy-all-read-expression-only]
                           [--skip-pretrim-fastqc]
                           [--skip-posttrim-fastqc]
                           [--skip-contamination-qc]
                           [--combine-contamination-qc]
                           [--contamination-qc-pairs CONTAMINATION_QC_PAIRS]
                           [--skip-alignment-qc] [--skip-umi-qc]
                           [--run-multiqc]
                           [-p PROJECT]

This script is used to generate input json files from the fastq_raw dir on gcp
for running rna-seq pipeline on GCP

optional arguments:
  -h, --help            show this help message and exit
  -g GCP_PATH, --gcp_path GCP_PATH
                        location of the submission batch directory in gcp that
                        contains the fastq_raw dir
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        output path, where you want the input jsons to be
                        written
  -r OUTPUT_REPORT_NAME, --output_report_name OUTPUT_REPORT_NAME
                        name of the output report to be written
  -u, --undetermined    Adding this flag will process undetermined FastQ files
                        if they exist. These are fastq files with prefix
                        "Undetermined_". If this flag isn't passed, items with
                        prefix "Undetermined_" will be removed
  -a {rat,human}, --organism {rat,human}
                        organism name, e.g. rat or human
  -v VERSION, --version VERSION
                        genome build or annotation release
  -n NUM_CHUNKS, --num_chunks NUM_CHUNKS
                        number of chunks to split the input files, should
                        always be <= number of input files
  -d DOCKER_REPO, --docker_repo DOCKER_REPO
                        Docker repository prefix containing the images used in
                        the workflow
  --release-manifest RELEASE_MANIFEST
                        complete release profile overriding reference and
                        Docker workflow inputs
  --runtime-profile RUNTIME_PROFILE
                        complete, explicitly selected CPU, RAM, and disk
                        override profile
  --star-disk-type {HDD,SSD}
                        STAR working-disk class; omit to retain the historical
                        HDD default
  -i, --index           add matched I1 FASTQs for UMI processing (enabled by
                        default)
  --legacy-all-read-expression-only
                        omit directional UMI molecule matrices and retain
                        historical all-read matrices only
  --skip-pretrim-fastqc
                        skip raw-read FastQC
  --skip-posttrim-fastqc
                        skip trimmed-read FastQC
  --skip-contamination-qc
                        skip globin, rRNA, and PhiX Bowtie2 screens
  --combine-contamination-qc
                        run full-depth globin, rRNA, and PhiX screens serially
                        on one worker
  --contamination-qc-pairs CONTAMINATION_QC_PAIRS
                        deterministically sample this many post-trim pairs for
                        the combined screens; 0 retains full depth
  --skip-alignment-qc
                        skip Picard duplicate/RNA metrics and chromosome
                        summaries
  --skip-umi-qc       skip directional UMI QC unless molecule expression
                        requires grouping
  --run-multiqc         emit legacy pre- and post-alignment MultiQC report
                        archives
  -p PROJECT, --project PROJECT
                        Project name on the google cloud platform

```

The metric-producing QC groups are independent and default to the historical
enabled behavior. Skipped metrics remain as empty fields in the stable QC
matrix. Cutadapt and STAR metrics remain available because those core
processing tasks always run. MultiQC is an opt-in reporting sidecar and does
not feed the QC matrix. The `--run-multiqc` compatibility mode requires both
FastQC groups and alignment QC.

The generator includes matched I1 FASTQs and directional UMI molecule-expression
matrices by default. The conventional all-read matrices are always retained.
Use `--legacy-all-read-expression-only` to omit only the additional molecule
matrices.

Default molecule-expression example:

```
python3 scripts/make_json_rnaseq.py -g gs://motrpac/rna-seq/test \
-o ./input_json \
-r rna-seq-test \
-a human \
-v gencode_v39 \
-n 1 \
-d us-docker.pkg.dev/motrpac-portal/rnaseq
```
