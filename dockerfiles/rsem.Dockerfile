FROM quay.io/biocontainers/rsem@sha256:05b9a0be934d1ac7fe894d965fbe45deff8d4fb1bbbfddb31b339541116b7420

# Keep the upstream stopping tolerance; allow slow EM fits more iterations.
ARG MAX_ROUND=20000
RUN curl -fsSL https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz -o /tmp/rsem.tar.gz \
    && echo '90e784dd9df8346caa2a7e3ad2ad07649608a51df1c69bfb6e16f45e611a40dc  /tmp/rsem.tar.gz' | sha256sum -c - \
    && tar -xzf /tmp/rsem.tar.gz -C /tmp \
    && cd /tmp/RSEM-1.3.3 \
    && test "$MAX_ROUND" -ge 20 && test "$MAX_ROUND" -le 100000 \
    && sed -i "s/^const int MAX_ROUND = 10000;$/const int MAX_ROUND = ${MAX_ROUND};/" EM.cpp \
    && sed -i 's/return (in>>/return bool(in>>/' SingleHit.h PairedEndHit.h \
    && x86_64-conda-linux-gnu-g++ -std=c++14 -O3 -ffast-math -DBOOST_NO_CXX98_FUNCTION_BASE=1 -I. -I/usr/local/include -c EM.cpp \
    && x86_64-conda-linux-gnu-g++ -std=c++14 -O3 -I. -I/usr/local/include -c SamHeader.cpp \
    && x86_64-conda-linux-gnu-g++ -pthread -Wl,-rpath,/usr/local/lib -L/usr/local/lib -o /usr/local/bin/rsem-run-em EM.o SamHeader.o -lhts -l:libz.so.1 \
    && mkdir -p /usr/local/share/rsem \
    && cp COPYING EM.cpp SingleHit.h PairedEndHit.h /usr/local/share/rsem/ \
    && rm -rf /tmp/rsem.tar.gz /tmp/RSEM-1.3.3
