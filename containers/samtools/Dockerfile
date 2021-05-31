FROM quay.io/comp-bio-aging/base:latest

ENV SAMTOOLS_VERSION 1.12

# update package manager & build essentials
RUN apt-get update \
    && apt-get install --yes build-essential wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev

# install samtools
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd / && rm -rf /tmp/samtools-${SAMTOOLS_VERSION}

CMD samtools