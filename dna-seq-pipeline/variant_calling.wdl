version development

workflow variant_calling {

    input {
        File bam
        File bai
        File referenceFasta
        File referenceFai
        Int threads = 16
        Int max_memory = 28
        String name
    }

      call smoove{
            input:
                bam = bam, bai = bai, reference = referenceFasta, reference_index = referenceFai, sample = name
        }

    call strelka2_germline{
        input:
            bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai, cores = threads
            #,indel_candidates = manta_germline_sv.manta_indel_candidates
    }


    output {
        File variants = strelka2_germline.variants
        File variantsIndex = strelka2_germline.variantsIndex
        File variants_SV = smoove.smooveVcf
        File results_SNP = strelka2_germline.results
        File results_SV = smoove.out
    }
}


task strelka2_germline{
    input {
        String runDir = "./strelka_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        #File? indel_candidates # ~{"--indelCandidates " + indel_candidates} \
        Boolean exome = false
        Boolean rna = false

        Int cores = 8
        Int max_memory = 28
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
        Int max_memory = 16
        Int max_cores = 8
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        smoove call \
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
        docker: "brentp/smoove@sha256:1bbf81b1c3c109e62c550783c2241acc1b10e2b161c79ee658e6abd011000c67"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{max_cores}"
    }
}