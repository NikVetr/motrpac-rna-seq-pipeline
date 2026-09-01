FROM ubuntu:20.04 as build

RUN apt-get update && \
    apt-get install -y --no-install-recommends gawk procps && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

COPY wdl/attach_umi/UMI_attach.awk /usr/local/src/UMI_attach.awk

RUN chmod 755 /usr/local/src/UMI_attach.awk
