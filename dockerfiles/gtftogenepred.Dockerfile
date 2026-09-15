FROM ubuntu:20.04@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a AS compiler

RUN apt-get update && \
    apt-get install -y --no-install-recommends wget ca-certificates && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

WORKDIR /usr/gtfToGenePred
RUN wget --progress=dot:giga https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/gtfToGenePred && \
    echo '306e7c1d8da3890bf813e4ae5c94a28bd9bd51ece9cb2b5ee2f51bf98e59f4fc  gtfToGenePred' | sha256sum -c -

FROM ubuntu:20.04@sha256:c664f8f86ed5a386b0a340d981b8f81714e21a8b9c73f658c4bea56aa179d54a AS build

RUN apt-get update && \
    apt-get install -y --no-install-recommends libkrb5-dev procps && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

COPY --from=compiler /usr/gtfToGenePred/gtfToGenePred /usr/local/bin/

RUN chmod 755 /usr/local/bin/gtfToGenePred
