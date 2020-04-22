FROM quay.io/comp-bio-aging/samtools:latest

ENV minimap2_version=2.17

RUN cd /tmp/ \
    && wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar -jxvf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && cd minimap2-${minimap2_version}_x64-linux \
    && cp minimap2 k8 paftools.js /usr/local/bin/

CMD minimap2