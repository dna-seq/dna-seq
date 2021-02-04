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
Before running the pipeline with a large genome (human or mouse) make sure you have 1-1.5 TB of free space. 

Prepare data
------------

[DVC](https://dvc.org/) is used for data management: it downloads annotations and can also be used to run some useful scripts.
To download the all the data and do some preprocessing use:
```bash
dvc repro
```
However, it may take quite a while as ensembl_vep_cache (which is required for VEP annotations) is >14GB. 
If you do not need it, you can also separately use dvc stages in dvc.yaml.
For instance, just to download reference genome, type:
```bash
dvc repro prepare_genome
```

Running services
----------------

To run [Cromwell server](https://cromwell.readthedocs.io/en/stable/), together with [cromwell-client](https://github.com/antonkulaga/cromwell-client) and mysql - run services with [Docker-Compose](https://docs.docker.com/compose/install/) 
If you do not have Docker installed, you can either install it yourself or use ubuntu_Script to in
For the convenience start.sh and stop bash scripts is put to the folder.
Cromwell services should already be working (I successfully tested with tests/hello-world.wdl on my Linux Mint machine)
\ To run the pipelines I recommend trying cromwell-client (deployed at port 8001 by default)
It can be run by:
```
docker compose up
```

In the docker-compose configuration the following configuration is assumed:
  
./data/databases/mysql folder for cromwell mysql

./data/cromwell-executions for cromwell execution cache

./data/cromwell-workflow-logs for cromwell execution logs

If you have another folder layout you have to change docker-compose.yml and config/cromwell/application.conf


Structure of the pipeline
-------------------------

The pipeline is in dna-seq-pipeline folder. 
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
    * /data/ensembl/102/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    * /data/ensembl/102/species/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
* OUTPUT folders created by the pipeline (no need to create them yourself)!: 
    * /data/gwas/anton/aligned - output of aligned data
    * /data/gwas/anton/variants - output for variants
    * /data/gwas/anton/vep - vep annotations
* ANNOTATION reference files (used only by vep_annotations.wdl):
    * /data/gwas/references/annotations - reference files for genetic variant annotations
    * /data/ensembl/102/plugins - git cloned+renamed https://github.com/Ensembl/VEP_plugins (note 102 - is Ensembl release number)
    * /data/ensembl/102/cache - folder to download Ensembl cache ( https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html ):
        VEP cache can be installed by:
        ```
        cd /data/ensembl/102/cache
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
* making haplausauras work and comparing with VEP
* eye and skin color annotations (stuff like http://mathgene.usc.es/snipper/eyeclassifier.html)
* FIX DUPLICATION OF COPY TASK
