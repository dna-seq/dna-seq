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

ENV SAMBABA_VERSION 0.8.0
WORKDIR /opt/sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v${SAMBABA_VERSION}/sambamba-${SAMBABA_VERSION}-linux-amd64-static.gz \
   && gunzip sambamba-${SAMBABA_VERSION}-linux-amd64-static.gz
RUN chmod +x sambamba-${SAMBABA_VERSION}-linux-amd64-static \
   && ln -s sambamba-${SAMBABA_VERSION}-linux-amd64-static sambamba
ENV PATH /opt/sambamba:$PATH
