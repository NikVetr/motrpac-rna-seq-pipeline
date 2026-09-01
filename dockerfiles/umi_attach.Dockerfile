FROM quay.io/biocontainers/gawk@sha256:14649372ffb4ac76e26ac7357262f2cdd3d731ace344118ad0ce79fc49f03a9d

COPY wdl/attach_umi/UMI_attach.awk /usr/local/src/UMI_attach.awk
