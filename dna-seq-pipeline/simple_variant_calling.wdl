version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow simple_variant_calling {

    input {
        File bam
        File bai
        String destination
        File referenceFasta
        File referenceFai
        Int threads# = 16
        Int max_memory# = 28
        String name
    }

    call smoove{
        input:
            bam = bam, 
            bai = bai, 
            reference = referenceFasta, 
            reference_index = referenceFai, 
            sample = name,
            max_memory = max_memory,
            max_cores = threads,
        }

    call strelka2_germline{
        input:
            bam = bam, 
            bai = bai, 
            reference_fasta = referenceFasta, 
            reference_fai= referenceFai, 
            cores = threads,
            max_memory = max_memory,
            #,indel_candidates = manta_germline_sv.manta_indel_candidates
    }

    call files.copy as copy_variants{
        input:
            destination = destination + "/variants",
            files =[
                   strelka2_germline.results,smoove.smooveVcf
                   ]
    }


    output {
        File variants = strelka2_germline.variants
        File variantsIndex = strelka2_germline.variantsIndex
        File results_SNP = copy_variants.out[0]
        File variants_SV = copy_variants.out[1]
        File results_SV = smoove.out
    }
}


task strelka2_germline{

    #deepvariant seems to be superior, but we keep strelka for the comparison purposes

    input {
        String runDir = "./strelka_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        #File? indel_candidates # ~{"--indelCandidates " + indel_candidates} \
        Boolean exome = false
        Boolean rna = false

        Int cores# = 8
        Int max_memory# = 28
    }

    command {
        configureStrelkaGermlineWorkflow.py \
        --bam ~{bam} \
        --ref ~{reference_fasta} \
        --runDir ~{runDir} \
        ~{true="--rna" false="" rna}
        ~{true="--exome" false="" exome} \
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{max_memory}
    }

    output {
        File variants = runDir + "/results/variants/variants.vcf.gz"
        File variantsIndex = runDir + "/results/variants/variants.vcf.gz.tbi"
        File results = runDir
    }

    runtime {
        docker: "quay.io/biocontainers/strelka:2.9.10--0"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{cores}"
    }
}


task smoove {
    input {
        File bam
        File bai
        File reference
        File reference_index
        String sample
        String outputDir = "./smoove"
        Int max_memory# = 16
        Int max_cores# = 8
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        smoove call \
        --genotype \
        --duphold \
        --outdir ~{outputDir} \
        --name ~{sample} \
        --fasta ~{reference} \
        ~{bam}
    }

    output {
        File out = outputDir
        File smooveVcf = outputDir + "/" + sample + "-smoove.vcf.gz"
    }

    runtime {
        docker: "brentp/smoove@sha256:d0d6977dcd636e8ed048ae21199674f625108be26d0d0acd39db4446a0bbdced"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{max_cores}"
    }
}
