FROM quay.io/comp-bio-aging/base:master

MAINTAINER Anton Kulaga

ENV PLATFORM=linux-x86_64

ENV BOWTIE2_VERSION=2.4.4

WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-${PLATFORM}.zip/download \
        -O bowtie2-${BOWTIE2_VERSION}-${PLATFORM}.zip && \
        unzip bowtie2-${BOWTIE2_VERSION}-${PLATFORM}.zip && \
        rm -rf bowtie2-${BOWTIE2_VERSION}-${PLATFORM}.zip && \
        ln -s bowtie2-${BOWTIE2_VERSION}-${PLATFORM} bowtie2

ENV PATH "$PATH:/opt/bowtie2"
CMD "/opt/bowtie2/bowtie2"