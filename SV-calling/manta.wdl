
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
