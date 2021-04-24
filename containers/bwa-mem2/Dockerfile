FROM quay.io/comp-bio-aging/base:latest
ENV VERSION=2.2.1

WORKDIR /opt
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-${VERSION}_x64-linux.tar.bz2 | tar jxf -
ENV PATH /opt/bwa-mem2-${VERSION}_x64-linux:$PATH

CMD /opt/bwa-mem2-${VERSION}_x64-linux/bwa-mem2