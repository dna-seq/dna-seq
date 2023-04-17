Just DNA-Seq pipelines
======================

Just a DNA-Seq analysis pipeline for personal and public genomes with some plugins and scripts on top of it.

Despite being expensive in the beginning, today genome sequencing is something accessible to all of us and costs roughly 400-800 dollars.
There are multiple sequencing and analysis proprietary services, however, their results are often based on proprietary databases and algorithms, and their predictions are often non-transparent.

Just DNA-Seq project was created primarily for transparency reasons: we wanted to understand what is happening. 
We also wanted to use the latest version of the tools as we discovered that, for example, DANTE-labs were using an outdated version of the genomes and GATK.
The project consists of multiple pipelines and scripts and can be either used separately or altogether.
Recently we started working on integration of [Plex Web 3 technology to run scientific workflows on decentralized infrastructure](https://github.com/labdao/plex).

About Just-DNA-Seq project, De-Sci and Web3
-------------------------------------------

This repository became a start of a larger [Just-DNA-Seq project](http://dna-seq.github.io/) with its own [Github Organization](https://github.com/dna-seq) and 23 source code repository.
Just-DNA-Seq project was funded by Longevity Gitcoin round in 2022 and reached most of its objective, including longevity and health risk genetic annotations, polygenic risk scores for longevity and multiple diseases, its scientific research paper is currently in a draft state.
While everybody could maintain their privacy by running both the pipelines and annotations on their own PC-s, the pipeline part requires some technical knowledge and a lot of RAM. 
Now we are working ( together with [LabDAO](https://labdao.xyz/) ) on adding PLEX support in our workflows to run them in a decentralized and privacy preserving ways.
We are also collaborating with [Genomes DAO](http://genomes.io) to add support of Just-DNA-Seq annotators and pipelines in DNA Vaults, owned and controlled by the user

Getting started
---------------

In the project we are using [WDL](https://openwdl.org/) (Workflow Description Language) pipelines as well as OakVar variant annotation system.
If you want to run the whole pipeline make sure you have at least 500GB or more of free space and >=16GB of RAM. 
All tools are dockerized, for this reason make sure that Docker is installed.

If genetic pipelines are something new for you, it can be useful to watch the Broad Institute [video introduction](https://www.youtube.com/watch?v=aTAQ2eA_iOc&feature=youtu.be&fbclid=IwAR0r2YeeJMEh2XFmat6OIEmbmGWXEvye3UYplvSheYFl7mJ1ijR65G0awLc) which explains WDL, Cromwell, and DNA-Seq pipelines.
Even though we do not use Broad-s GATK pipeline and mix our tools a bit differently (for example we use [DeepVariant](https://academic.oup.com/bioinformatics/article/36/24/5582/6064144) for variant calling), the video explains some useful concepts in genomic analysis.
For the users with only high-school knowledge of biology, I would also recommend taking any free biology 101 or genetics 101 course ( https://www.edx.org/course/introduction-to-biology-the-secret-of-life-3 is a good example)

For gene annotations we use [OakVar](https://rkimoakbioinformatics.github.io/oakvar/) as well as VEP (as an alternative solution).
OakVar is included in the conda environment. 


Install conda environment
-------------------------
Annotation modules and dvc are included in the conda environment that goes with the project.
The environment can be setup either with Anaconda or micromamba (superior version of conda).
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
micromamba activate dna-seq
```

The instructions above are provided for Linux and MacOS (note: in MacOS you have to install wget). 
For Windows you can either install Linux Subsystem or use Windows version of anaconda.

Prepare data
------------

[DVC](https://dvc.org/) is used for data management: it downloads annotations and can also be used to run some useful scripts.
DVC is included to the gwas conda environment described in environment.yaml file.

In dvc.yaml there are tasks required to setup the project. For instance, to download reference genome, type:
```bash
dvc repro prepare_genome
```
Of course, you can try to download all the data with:
```bash
dvc repro
```
However, it may take quite a while as ensembl_vep_cache (which is required for VEP annotations) is >14GB. 
And it may happen that OakVar will be enough for your needs.
In the Future, we plan to focus on OakVar leaving VEP as a legacy annotation system.


Running services
----------------

To run [Cromwell server](https://cromwell.readthedocs.io/en/stable/), together with [cromwell-client](https://github.com/antonkulaga/cromwell-client) and mysql - run services with [Docker-Compose](https://docs.docker.com/compose/install/) 
If you do not have Docker installed, you can either install it yourself or use ubuntu_Script in the bin folder.
To run the pipelines I recommend trying cromwell-client (deployed at port 8001 by default).
It can be run by:
```bash
docker compose up
```

In the docker-compose configuration the following configuration is assumed:
  
./data/databases/mysql folder for cromwell mysql

./data/cromwell-executions for cromwell execution cache

./data/cromwell-workflow-logs for cromwell execution logs

If you have another folder layout you have to change docker-compose.yml and config/cromwell/application.conf.


Structure of the pipeline
-------------------------

The pipeline is in dna-seq-pipeline folder. It also actively uses wdl tasks and subpipelines from https://github.com/antonkulaga/bioworkflows.
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
    * /data/ensembl/103/cache - folder to download Ensembl cache ( https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html )
        
Genomes to play with
------------------

I tested the pipeline on my personal genome sequenced by Dante. 
If you do not have your own genome at disposal, you can try any of the public ones, for example, you can download WGS fastq files from https://www.personalgenomes.org.uk/data/.
For a quick test of all tools consider having small test fastq-s (example of such test is in test.json).
If you have a bam file with input we provide bam_to_fastq pipeline to extract fastq-s from it.


Different options of running the pipelines
-----------------

There are three major ways of running any of the pipelines in the repository:
* with [CromwellClient](https://github.com/antonkulaga/cromwell-client) and Cromwell in server mode (recommended). Note: when using pipelines with multiple wdl files, do not forget to upload subworkflow files as dependencies.
* directly from Swagger API with Cromwell in a server mode: similar to running with CromwellClient but instead of the Client swagger server API is used.
* with Cromwell or any other WDL-compartible tool in the console. Documented at [Official cromwell documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/#step-3-running-the-workflow).

Genome annotations
==================

There are two alternative annotation tools: VEP and [OakVar](https://rkimoakbioinformatics.github.io/oakvar/). VEP is more established and oldfashioned. OakVar is newer and more user-friendly. I recommend starting from OakVar.

OakVar annotations
======================

OakVar is included to the environment.
Before starting, it is recommended to install annotation modules of your interest.
There is a dvc stage for the default modules:
```bash
dvc repro install_oakvar
```

VEP annotations
---------------

The most important part of the pipeline is VEP annotations.
To make it work you must download Ensembl cache and Ensembl plugins.
Right now VEP annotations require additional files, for this reason, they are provided as a separate pipeline.
Currently DVC resolves most of the files and does additional preprocessing with:
```
dvc repro prepare
```

Development plan
----------------

The pipeline still requires some technical skills to run, we plan to improve ease of use and stream-line it a bit.
One of the most important parts is variant filtering and annotations. 
OakVar does a great job of installing a huge number of annotation sources. 
However, its output requires some biological skills to read, we are working now on:
* reporting plugins to make reports for pre-selected genes
* longevity plugin for analysis of gene variants associated with longevity.
* mitochondrial variant calling
* longevity annotations

Links
----------------

* Documentation:  https://just-dna-seq.readthedocs.io/en/oakvar/
* Discord: https://discord.gg/5WU6aSANXy
