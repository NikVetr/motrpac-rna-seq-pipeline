# Scripts

### `make_json_rnaseq.py`

Generates the input configuration file required to run the rna-seq pipeline.

- Requires Python `>3.6.9`
- Install required packages by running `pip3 install -r scripts/requirements.txt`

```
usage: make_json_rnaseq.py [-h] -g GCP_PATH -o OUTPUT_PATH
                           -r OUTPUT_REPORT_NAME [-u] -a {rat,human}
                           -v {rn6,rn7,rn8,gencode_v39,gencode_v47}
                           -n NUM_CHUNKS
                           [--sample-list SAMPLE_LIST |
                            --exclude-sample-list EXCLUDE_SAMPLE_LIST]
                           [-d DOCKER_REPO]
                           [--release-manifest RELEASE_MANIFEST]
                           [--runtime-profile RUNTIME_PROFILE]
                           [--star-disk-type {HDD,SSD}] [-i]
                           [--legacy-all-read-expression-only]
                           [--retain-all-read-expression]
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
  --sample-list SAMPLE_LIST
                        process only the exact sample prefixes listed one per
                        line
  --exclude-sample-list EXCLUDE_SAMPLE_LIST
                        exclude the exact sample prefixes listed one per line
  -d DOCKER_REPO, --docker_repo DOCKER_REPO
                        Docker repository prefix containing the images used in
                        the workflow
  --release-manifest RELEASE_MANIFEST
                        complete release profile overriding reference and
                        Docker workflow inputs
  --runtime-profile RUNTIME_PROFILE
                        complete CPU, RAM, and disk-floor override profile
  --star-disk-type {HDD,SSD}
                        STAR working-disk class; omit to retain the historical
                        HDD default
  -i, --index           add matched I1 FASTQs for UMI processing (enabled by
                        default)
  --legacy-all-read-expression-only
                        omit directional UMI molecule matrices and retain
                        historical all-read matrices only
  --retain-all-read-expression
                        also run and emit non-UMI-deduplicated matrices under
                        all_read names
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
matrices by default. They use the canonical expression filenames and are the
only expression branch run by default. Use `--retain-all-read-expression` to
also emit the non-UMI-deduplicated matrices under `all_read_*` names, or use
`--legacy-all-read-expression-only` to make the historical all-read matrices
canonical instead. The two switches cannot be combined.

Runtime profiles set fixed CPU and memory requests plus a minimum STAR scratch
size. After Cutadapt, the workflow uses each sample's exact surviving pair
count to raise STAR scratch independently when needed, so one generated JSON
can contain samples in different disk tiers.

For a bounded pilot, write the exact sample prefixes (the FASTQ basename before
`_R1.fastq.gz`) one per line and pass `--sample-list`. Reuse that same manifest
with `--exclude-sample-list` when generating the remaining samples. Selection
is within the single GCS prefix passed to `--gcp_path`; one JSON may contain all
selected samples when `--num_chunks 1` is used.

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
