FROM quay.io/biocontainers/python@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2

COPY wdl/collect_qc_metrics/rnaseq_qc.py /usr/local/src/rnaseq_qc.py
