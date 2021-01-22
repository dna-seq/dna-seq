#!/bin/bash

#env
ens_ver=102
ens_url="ftp://ftp.ensembl.org/pub/release-${ens_ver}"
ens_dir=/data/gwas/homo_sapiens_ensembl_${ens_ver}
GRCh38="Homo_sapiens.GRCh38"
GRCh38_primary_asm="${GRCh38}.dna.primary_assembly"

echo cromwell/db folders...
mkdir -p /data/databases/mysql
mkdir -p /data/cromwell-workflow-logs
mkdir -p /data/cromwell-executions

#REFERENCE genome

mkdir -p ${ens_dir}
	#Reference Genome
if ! [ -f "${ens_dir}/${GRCh38_primary_asm}.fa" ]; then
	echo Building Ensembl Reference Genome from ${ens_url}/fasta/homo_sapiens/dna/${GRCh38_primary_asm}.fa.gz
	wget -qO- ${ens_url}/fasta/homo_sapiens/dna/${GRCh38_primary_asm}.fa.gz \
	| gunzip -c >  ${ens_dir}/${GRCh38_primary_asm}.fa
fi
	#Index of Reference genome
#if ! [ -f "${ens_dir}/${GRCh38_primary_asm}.fa.fai" ]; then
#	echo Building FAI Index of Reference genome from ${ens_dir}/${GRCh38_primary_asm}.fa
#	bwa index ${ens_dir}/${GRCh38_primary_asm}.fa; 
#	samtools faidx ${ens_dir}/${GRCh38_primary_asm}.fa 
#fi

	#Dictionary of Reference genome
#if ! [ -f "${ens_dir}/${GRCh38_primary_asm}.dict" ]; then
#	echo Building Dictionary of Reference genome from ${ens_dir}/${GRCh38_primary_asm}.fa
#	java -jar /usr/local/bin/picard.jar CreateSequenceDictionary -R ${ens_dir}/${GRCh38_primary_asm}.fa -O ${ens_dir}/${GRCh38_primary_asm}.dict 
#fi

	#Gene Annotation
#if ! [ -f "${ens_dir}/${GRCh38}.${ens_ver}.gtf" ]; then
#	echo Building GTF Gene Annotation of Reference genome from ${ens_url}/gtf/homo_sapiens/${GRCh38}.${ens_ver}.gtf.gz
#	wget -qO- ${ens_url}/gtf/homo_sapiens/${GRCh38}.${ens_ver}.gtf.gz | gunzip -c > ${ens_dir}/${GRCh38}.${ens_ver}.gtf 
#fi

#Ensembl cache
if ! [ -d "/data/ensembl/${ens_ver}/cache/homo_sapiens" ]; then
	echo Building Ensembl cachce from ${ens_url}/variation/vep/homo_sapiens_vep_${ens_ver}_GRCh38.tar.gz\n
	mkdir -p /data/ensembl/${ens_ver}/cache
	#indexed
	wget -qO- ${ens_url}/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz | tar zxvf -
	#non-indexed
	#wget -qO- ${ens_url}/variation/vep/homo_sapiens_vep_${ens_ver}_GRCh38.tar.gz | tar zxvf -
	
fi

#Ensembl plugins 
if ! [ -d "/data/ensembl/${ens_ver}/plugins" ]; then
	echo Building Ensembl plugins from https://github.com/Ensembl/VEP_plugins\n
	mkdir -p /data/ensembl/${ens_ver}/plugins
	git clone https://github.com/Ensembl/VEP_plugins /data/ensembl/${ens_ver}/plugins
fi

echo Complete!
