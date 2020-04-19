version development

workflow Lumpy {
    input {
        File ref
        File ref_fai
        File bam
        File bai
        String name
        String destination

    }

    call Smoove{
        input:
            bam = bam, bai = bai, reference = ref, reference_index = ref_fai, sample = name
    }

    call copy {
        input: destination = destination + "/lumpy", files = [Smoove.out]
    }

    output {
        File out = copy.out[0]
    }
}


task Smoove {
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

task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
          do
              value=$(basename ~{"$"}i)
              echo ~{where}/~{"$"}value
          done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}