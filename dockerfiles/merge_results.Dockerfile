FROM quay.io/biocontainers/python@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2

COPY wdl/merge_results/consolidate_qc_report.py /usr/local/src/consolidate_qc_report.py
COPY wdl/merge_results/merge_fc.py /usr/local/src/merge_fc.py
COPY wdl/merge_results/merge_rsem.py /usr/local/src/merge_rsem.py
