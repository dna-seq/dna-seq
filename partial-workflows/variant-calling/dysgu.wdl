version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow dysgu {
    input {
        File ref
        File ref_fai
        File bam
        File bai
        String name
        String destination

    }

    call dysgu {
        input:
            bam = bam, bai = bai, reference = ref, reference_index = ref_fai, sample = name
    }

    call files.copy as copy {
        input: destination = destination, files = [dysgu.out, dysgu.out_bam]
    }

    output {
        File out = copy.out[0]
    }
}


task dysgu {
    input {
        File bam
        File bai
        File reference
        File reference_index
        String sample
        String outputDir = "./dysgu"
        Int max_memory = 16
        Int max_cores = 8
    }

    String bam_input = "inputs/" + basename(bam)
    String bai_input = "inputs/" + basename(bai)
    String reference_input = "inputs/" + basename(reference)
    String reference_index_input = "inputs/" + basename(reference_index)

    # --remap False --metrics True

    command {
        set -e
        mkdir -p inputs
        ln -s ~{bam} ~{bam_input}
        ln -s ~{bai} ~{bai_input}
        ln -s ~{reference} ~{reference_input}
        ln -s ~{reference_index} ~{reference_index_input}
        dysgu run --diploid False --mode pe ~{reference_input} ~{outputDir} ~{bam_input} > ~{sample}.vcf
    }

    output {
        File out = sample + ".vcf"
        File out_bam = outputDir + "/" + sub(basename(bam), ".bam", ".dysgu_reads.bam")
    }

    runtime {
        docker: "kcleal/dysgu:latest"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{max_cores}"
    }
}