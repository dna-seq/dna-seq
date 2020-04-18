version development

workflow variant_calling {

    input {
        File bam
        File bai
        File referenceFasta
        File referenceFai
        Int threads = 16
        Int max_memory = 28
    }

    call manta_germline_sv{
            input:
                bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai, cores = threads
        }

    call strelka2_germline{
        input:
            bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai, cores = threads,
            indel_candidates = manta_germline_sv.manta_indel_candidates
    }

    call vep_annotation{
        input: vcf = strelka2_germline.variants
    }

    call vep_annotation as manta_vep_annotation{
        input: vcf = manta_germline_sv.manta_SV
    }

    output {
        File variants = strelka2_germline.variants
        File variantsIndex = strelka2_germline.variantsIndex
        File results_SNP = strelka2_germline.results
        File manta_SV = manta_germline_sv.manta_SV
        File manta_SV_index = manta_germline_sv.manta_SV_index

        File annotations = vep_annotation.out
        File vep_summary = vep_annotation.summary

        File manta_annotations = manta_vep_annotation.out
        File manta_vep_summary = manta_vep_annotation.summary
    }
}


task manta_germline_sv {
    input {
        String runDir = "./manta_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        Boolean exome = false
        Boolean rna = false

        Int cores = 8
        Int max_memory = 28
    }

    command {
        set -e
        configManta.py \
        ~{"--normalBam " + bam} \
        --referenceFasta ~{reference_fasta} \
        --runDir ~{runDir} \
        ~{true="--exome" false="" exome}
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{max_memory}
    }

    output {
        File manta_SV = runDir + "/results/variants/diploidSV.vcf.gz"
        File manta_SV_index = runDir + "/results/variants/diploidSV.vcf.gz.tbi"
        File? manta_indel_candidates = runDir + "/results/variants/candidateSmallIndels.vcf.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/manta@sha256:81a24337616d9bceaf5afdeed1c5b70ba5f238194389717ad7ab2f188de17bfa"#:1.6.0--py27_0"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{cores}"
    }
}

task strelka2_germline{
    input {
        String runDir = "./strelka_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        File? indel_candidates
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
        ~{"--indelCandidates " + indel_candidates} \
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



task vep_annotation {
    input {
        File vcf
        String name = "variant_effect_output.tsv"
        String species = "human"
        Int threads = 8
        Boolean offline = true
    }

    #TODO add VEP plugins http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html

    command {
        vep \
        --input_file ~{vcf} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything ~{if(offline) then "--offline" else ""} \
        --gene_phenotype --biotype --uniprot --symbol --allele_number --total_length --allele_number --regulatory --af
    }

    runtime {
        docker: "ensemblorg/ensembl-vep:release_99.2"
    }

    output {
        File out = name
        File summary = "variant_effect_output.txt_summary.html"
    }
}