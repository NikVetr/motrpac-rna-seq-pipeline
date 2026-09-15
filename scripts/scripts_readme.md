# Scripts

## Input generation

`make_json_rnaseq.py` generates workflow input JSONs from paired FASTQs under
one GCS prefix. Use Python >= 3.10 and install `scripts/requirements.txt`.
Run `python3 scripts/make_json_rnaseq.py --help` for the complete argument list.

```bash
mkdir -p input_json
python3 scripts/make_json_rnaseq.py \
  -g gs://BUCKET/SUBMISSION -o input_json -r cohort \
  -a human -v gencode_v50 -n 1 \
  --release-manifest config/release-profiles/human-gencode-v50.json \
  --runtime-profile config/backends/gcp/runtime-human-v50-full-candidate-v1.json
```

Choose a matched release and runtime profile for the organism and annotation.
See [`config/release-profiles/`](../config/release-profiles/) and
[`config/backends/gcp/`](../config/backends/gcp/) for available profiles.
Forward-stranded molecule expression is canonical when matched I1 reads are
present. `--allow-missing-umis` permits samples without I1 to use all-read
expression, recording the skip reason and `not_deduplicated` covariate.
`--all-read-expression-only` selects all-read counting for the whole submission;
`--retain-all-read-expression` adds a secondary all-read pass. These last two
options cannot be combined, and neither changes strandedness.

Use `--sample-list` or `--exclude-sample-list` with exact sample prefixes, one
per line, to select or exclude samples under the supplied GCS prefix. With
`--num_chunks 1`, all selected samples go into one workflow input JSON.

## QC and resources

FastQC, contamination, alignment and UMI QC groups are independently selectable.
Skipped metrics remain empty in the QC matrix; Cutadapt and STAR metrics are
always available. MultiQC archives are opt-in via `--run-multiqc` and require
both FastQC groups and alignment QC. `--combine-contamination-qc` shares one
worker across the three screens; `--contamination-qc-pairs` optionally limits
them to a deterministic post-trim pair sample.

Runtime profiles set CPU requests and RAM/disk floors. The workflow raises
STAR scratch from post-trim pair counts, UMI scratch from STAR BAM sizes, and
RSEM RAM/scratch from the transcriptome BAM entering quantification. See
[cohort provisioning](../docs/cohort-provisioning.md) for rules, deployment,
metadata preparation and retention.

## Maintenance

Build and optionally publish images with `build_dockerfiles.sh`; see the
[container instructions](../README.md#building-and-updating-containers).
`validate_jsons.py first.json second.json` compares JSON values for equality.
Use WOMtool to validate workflow inputs, as described in
[local testing](../README.md#validating-inputs-and-running-tests).
