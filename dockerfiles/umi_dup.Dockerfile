FROM quay.io/biocontainers/umi_tools@sha256:94c7cd9a713157affe93d3f1fa60e60d35a6385adc6b419d5f73c68eea8a54e8

COPY wdl/umi_dup/prepare_umi_bam.py /usr/local/src/prepare_umi_bam.py
COPY wdl/umi_dup/summarize_umi_tools.py /usr/local/src/summarize_umi_tools.py
COPY wdl/umi_dup/propagate_molecule_qnames.py /usr/local/src/propagate_molecule_qnames.py
COPY wdl/umi_dup/summarize_molecule_expression.py /usr/local/src/summarize_molecule_expression.py

RUN chmod 0555 \
      /usr/local/src/prepare_umi_bam.py \
      /usr/local/src/summarize_umi_tools.py \
      /usr/local/src/propagate_molecule_qnames.py \
      /usr/local/src/summarize_molecule_expression.py
