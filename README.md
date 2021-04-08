Simple DNA-Seq
==============
A simple pipeline to assemble personal and public genomes

Description
-----------

The pipeline is based on a WDL (Workflow Description Language) standard.
Broad Institute provides a nice [video introduction](https://www.youtube.com/watch?v=aTAQ2eA_iOc&feature=youtu.be&fbclid=IwAR0r2YeeJMEh2XFmat6OIEmbmGWXEvye3UYplvSheYFl7mJ1ijR65G0awLc) which explains WDL, Cromwell and DNA-Seq pipelines. 
For users with only high-school knowledge of biology I would also recommend taking any free biology 101 or genetics 101 course ( https://www.edx.org/course/introduction-to-biology-the-secret-of-life-3 is a good example)

We do not use Broad-s GATK pipeline (because we use DeepVariant as a variant caller) but common tools are similar. 
All tools are dockerized, for this reason make sure that docker is installed. 
Before running the pipeline with large genomes (human or mouse) make sure you have around 1 TB or more of free space.

Install conda environment
-------------------------
Annotation modules and dvc are included in the conda environment that goes with the project.
The environment can be setup either with anaconda or micromamba (superior version of conda).
Here is installation instruction for micromamba:
```
wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
```
We can use ./micromamba shell init ... to initialize a shell (.bashrc) and a new root environment in ~/micromamba:
```
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
```
To create a micromamba environment use:
```
micromamba create -f environment.yaml
micromamba activate gwas
```

Prepare data
------------

[DVC](https://dvc.org/) is used for data management: it downloads annotations and can also be used to run some useful scripts.
DVC is included to the gwas conda environment described in environment.yaml file

In dvc.yaml there are tasks required to setup the project. For instance, to download reference genome, type:
```bash
dvc repro prepare_genome
```
Of course, you can try to download all the data with:
```bash
dvc repro
```
However, it may take quite a while as ensembl_vep_cache (which is required for VEP annotations) is >14GB. And it may happen that OpenCravat will be enough for you needs.


Running services
----------------

To run [Cromwell server](https://cromwell.readthedocs.io/en/stable/), together with [cromwell-client](https://github.com/antonkulaga/cromwell-client) and mysql - run services with [Docker-Compose](https://docs.docker.com/compose/install/) 
If you do not have Docker installed, you can either install it yourself or use ubuntu_Script in the bin folder.
To run the pipelines I recommend trying cromwell-client (deployed at port 8001 by default)
It can be run by:
```bash
docker compose up
```

In the docker-compose configuration the following configuration is assumed:
  
./data/databases/mysql folder for cromwell mysql

./data/cromwell-executions for cromwell execution cache

./data/cromwell-workflow-logs for cromwell execution logs

If you have another folder layout you have to change docker-compose.yml and config/cromwell/application.conf


Structure of the pipeline
-------------------------

The pipeline is in dna-seq-pipeline folder. It also actively uses wdl tasks and subpipelines from https://github.com/antonkulaga/bioworkflows
There dna_seq_pipeline.wdl is the main workflow, all others should be provided as dependencies.
The pipeline uses Deepvariant, Strelka2 and Smoove for variant calling and VEP for variant annotations.
It is also possible to run dependencies as separate workflows.

Example JSONs
------------
Example json inputs are provided with the parameters that I used to process my own genome.
In the examples I had the following structure.
This structure is the same as in ./data subfolder of this repo, feel free to modify it according to locations of your files:
* /data/gwas/anton - person's folder with:
* /data/gwas/anton/fastq - INPUT folder
* REFERENCE genome (downloaded from latest Ensembl release):
    * /data/ensembl/103/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    * /data/ensembl/103/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
* OUTPUT folders created by the pipeline (no need to create them yourself, when you run the pipeline it will create folders for the output)!: 
    * /data/gwas/anton/aligned - output of aligned data
    * /data/gwas/anton/variants - output for variants
    * /data/gwas/anton/vep - vep annotations
* ANNOTATION reference files (used only by vep_annotations.wdl):
    * /data/gwas/references/annotations - reference files for genetic variant annotations
    * /data/ensembl/103/plugins - git cloned+renamed https://github.com/Ensembl/VEP_plugins (note 102 - is Ensembl release number)
    * /data/ensembl/103/cache - folder to download Ensembl cache ( https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html ):
        VEP cache can be installed by:
        ```
        cd /data/ensembl/103/cache
        curl -O ftp://ftp.ensembl.org/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz
        tar xzf homo_sapiens_vep_102_GRCh38.tar.gz
        cd $HOME/.vep       
        ```

Genomes to play with
------------------

I tested the pipeline on my personal genome sequenced by Dante. 
If you do not have your own genome in disposal, you can try any of the public ones, for example you can download WGS fastq files from https://www.personalgenomes.org.uk/data/
For quick test of all tools consider having small test fastq-s (example of such test is in test.json)
If you have a bam file with input we provide bam_to_fastq pipeline to extract fastq-s from it.


Different options of running the pipelines
-----------------

There are three major ways of running any of the pipelines in the repository:
* with [CromwellClient](https://github.com/antonkulaga/cromwell-client) and Cromwell in server mode (recommended). Note: when using pipelines with multiple wdl files, do not forget to upload subworkflow files as dependencies.
* directly from Swagger API with Cromwell in a server mode: similar to running with CromwellClient but instead of the Client swagger server API is used.
* with Cromwell or any other WDL-compartible tool in the console. Documented at [Official cromwell documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/#step-3-running-the-workflow)

Genome annotations
==================

There are two alternative annotation tools: VEP and [OpenCravat](https://opencravat.org/). VEP is more established and oldfasioned, opencravat is newer and more user-friendly I recommend to start from opencravat.

Opencravat annotations
======================

Opencravat is included to the environment.
Before starting, it is recommended to install annotation modules of your interested.
There is a dvc stage for the default modules:
```bash
dvc repro install_opencravat
```

VEP annotations
---------------

The most important part of the pipeline is VEP annotations.
To make it work you must download Ensembl cache and Ensembl plugins.
Right now VEP annotations require additional files, for this reason they are provided as separate pipeline.
Curently DVC resolves most of the files and does additional preprocessing with:
```
dvc repro prepare
```

TODOs
-----

This repository is a work in progress, so I list todos:
* mitochondrial variant calling
* longevity annotations
* docs on how to configure VEP
* docs on opencravat  
* eye and skin color annotations (stuff like http://mathgene.usc.es/snipper/eyeclassifier.html )