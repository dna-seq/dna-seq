Simple DNA-Seq
==============
A simple pipeline to assembl personal and public genomes

Running services
----------------

The pipeline is based on WDL standard.
Broad Institute provides a nice [video introduction](https://www.youtube.com/watch?v=aTAQ2eA_iOc&feature=youtu.be&fbclid=IwAR0r2YeeJMEh2XFmat6OIEmbmGWXEvye3UYplvSheYFl7mJ1ijR65G0awLc) which explains WDL, Cromwell and DNA-Seq pipelines.
All tools are dockerized. 
To run [Cromwell server](https://cromwell.readthedocs.io/en/stable/), together with [cromwell-client](https://github.com/antonkulaga/cromwell-client) and mysql - run services with [Docker-Compose](https://docs.docker.com/compose/install/) 
For the convenience start.sh and stop bash scripts is put to services folder.

Cromwell services should already be working (I successfully tested with tests/hello-world.wdl on my Linux Mint machine)
To run the pipelines I recommend trying my cromwell-client (deployed at port 8001 by default)
In the docker-compose configuration the following configuration is assumed:
  
/data/databases/mysql folder for cromwell mysql

/data/cromwell-executions for cromwell execution cache

/data/cromwell-workflow-logs for cromwell execution logs

If you have another folder layout you have to change docker-compose.yml and config/cromwell/application.conf


Structure of the pipeline
-------------------------

The pipeline is in dna-seq-pipeline folder. There dna_seq_pipeline.wdl is the main workflow, all others should be provided as dependencies


Different options of running the pipelines
-----------------

There are three major ways of running any of the pipelines in the repository:
* with [CromwellClient](https://github.com/antonkulaga/cromwell-client) and Cromwell in server mode (recommended)
* directly from Swagger API with Cromwell in a server mode: similar to running with CromwellClient but instead of the Client swagger server API is used.
* with Cromwell or any other WDL-compartible tool in the console. Documented at [Official cromwell documentation](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/#step-3-running-the-workflow)

VEP annotations
---------------

The most important part of the pipeline is VEP annotations.
To make it work you must download Ensembl cache and Ensembl plugins.

TODOs
-----

This repository is a work in progress, so I list todos:
* longevity annotations
* docs on how to configure VEP (maybe dvc-based resolved cache)
* making haplausauras work and comparing with VEP
* eye and skin color annotations (stuff like http://mathgene.usc.es/snipper/eyeclassifier.html)
* moving to DeepVariant as a main variant caller