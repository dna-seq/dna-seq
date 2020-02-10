Simple DNA-Seq
==============
A simple pipeline for personal use to analyze personal and public data

Running services
----------------

The pipeline is based on WDL standard.
Broad Institute provides a nice [video introduction](https://www.youtube.com/watch?v=aTAQ2eA_iOc&feature=youtu.be&fbclid=IwAR0r2YeeJMEh2XFmat6OIEmbmGWXEvye3UYplvSheYFl7mJ1ijR65G0awLc) which explains WDL, Cromwell and DNA-Seq pipelines.
All tools are dockerized. 
To run [Cromwell server](https://cromwell.readthedocs.io/en/stable/), together with [cromwell-client](https://github.com/antonkulaga/cromwell-client) and mysql - run services with [Docker-Compose](https://docs.docker.com/compose/install/) 
For the convenience start.sh and stop bash scripts is put to services folder. 
Cromwell services should already be working (I successfully tested with tests/hello-world.wdl on my Linux Mint machine)
To run the pipelines I recommend to try my cromwell-client (deployed at port 8001 by default)

Structure of the pipeline
-------------------------
dna_seq_pipeline.wdl is the main workflow, all others should be provided as dependencies
