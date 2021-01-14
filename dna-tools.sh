#!/bin/bash
if ! [ $(id -u) = 0 ]; then
   echo "Must be root!"
   exit 1
fi
apt-get update
#transport
apt-get install rsync tree
#code 
apt-get install gdb openjdk-8-jre-headless openjdk-8-jdk gnutls-bin
#libs
apt-get install libncurses-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libidn11-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev libssl-dev libmpfr-dev
#genetic
apt-get install bwa samtools bcftools

if ! [ -f "/usr/local/bin/samtools" ]; then
	mkdir -p /data/source/samtools
	cd /data/source/samtools
	wget -qO- https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 | tar jxvf - --strip-components 1
	./configure --prefix=/usr/local/
	make && make install
	cd -
fi
samtools --version
if ! [ -f "/usr/local/bin/bcftools" ]; then
	mkdir -p /data/source/bcftools
	cd /data/source/bcftools
	wget -qO- https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 | tar jxvf - --strip-components 1
	./configure --prefix=/usr/local/
	make && make install
	cd -
fi
bcftools --version
if ! [ -f "/usr/local/bin/bam2fastq" ]; then
	mkdir -p /data/source/bam2fastq
	git clone --recursive https://github.com/jts/bam2fastq /data/source/bam2fastq
	cd /data/source/bam2fastq
	make && mv ./bam2fastq /usr/local/bin/bam2fastq
	cd -
fi
bam2fastq --version
if ! [ -f "/usr/local/bin/picard.jar" ]; then
	mkdir -p /data/source/picard 
	git clone https://github.com/broadinstitute/picard.git /data/source/picard
	cd /data/source/picard
	./gradlew shadowJar && mv /data/source/picard/build/libs/picard.jar /usr/local/bin/picard.jar
	cd -
fi
java -jar /usr/local/bin/picard.jar CreateSequenceDictionary --version
