version development

workflow alignment {
    input {
        Array[File]+ reads
        File reference
        String name
        Int align_threads = 12
        Int sort_threads = 12
        Int max_memory_gb = 36
        Int coverage_sampling = 1000
        String? gencore_quality
    }

    call minimap2 {
        input:
            reads = reads,
            reference = reference,
            name = name,
            threads = align_threads,
            max_memory = max_memory_gb
    }

    call sambamba_sort {
        input:
            bam = minimap2.bam,
            threads = sort_threads
    }

    call gencore{
        input:
            reference = reference,
            sorted_bam = sambamba_sort.out,
            name = name,
            quality  = gencore_quality,
            coverage_sampling = coverage_sampling,
            max_memory = max_memory_gb
    }


    output {
       File bam = gencore.bam
       File bai = gencore.bai
       File html = gencore.html
       File json = gencore.json
    }
}


task minimap2 {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
        Int max_memory
    }

    command {
        minimap2 -ax sr  -t ~{threads} -2 ~{reference} ~{sep=' ' reads} | samtools view -bS - > ~{name}.bam
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker: "quay.io/comp-bio-aging/minimap2@sha256:f5d43a4d857fc56bfa4e98df1049a7b9c8af0f1bf604580eb074953a00b455cd" #latest
        maxRetries: 2
      }

    output {
      File bam = name + ".bam"
    }
}

task sambamba_sort{
    input {
        File bam
        Int threads
        Int gb_per_thread = 3
    }

    String name = basename(bam, ".bam")

    command {
       ln -s ~{bam} ~{basename(bam)}
       sambamba sort -m ~{gb_per_thread}G -t ~{threads} -p ~{basename(bam)}
       mv -f ~{name}.sorted.bam.bai ~{name}.sorted.bai
    }

    runtime {
        docker: "quay.io/biocontainers/sambamba@sha256:8aa120d440ff188d447eaa0e6d5cac82bd9e35bfa42d5c7857c401736629c299" #:0.7.1--h148d290_2
        maxRetries: 1
        docker_memory: "~{gb_per_thread * (threads+1)}G"
        docker_cpu: "~{threads+1}"
        docker_swap: "~{gb_per_thread * (threads+1) * 2}G"
      }

    output {
      File out = name + ".sorted.bam"
      File bai = name + ".sorted.bai"
    }
}


task gencore {
    input {
        File reference
        File sorted_bam
        String name
        Int max_memory
        Int supporting_reads = 1
        Float ratio_threshold = 0.8
        String? quality #"--high_qual"
        Int coverage_sampling = 1000
    }
    command {
        gencore --coverage_sampling ~{1000} --ratio_threshold=~{ratio_threshold} -s ~{supporting_reads} ~{quality} -i ~{sorted_bam} -o ~{name}.bam -r ~{reference}
        samtools index ~{name}.bam  ~{name}.bai
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker: "quay.io/comp-bio-aging/gencore@sha256:8bd1ef4984300abf194d94e1885c6642a72803fb8f81439205077c8db9f31103"
    }

    output {
        File bam = name + ".bam"
        File bai = name + ".bai"
        File html = "gencore.html"
        File json = "gencore.json"
    }
}



task coverage {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
        bedtools genomecov -bg -ibam ~{bam} > ~{name}.bedgraph
    }

     runtime {
            docker: "quay.io/biocontainers/bedtools@sha256:02e198f8f61329f9eafd1b9fc55828a31020b383403adec22079592b7d868006" #2.29.2--hc088bd4_0
            maxRetries: 2
          }

    output {
        File out = name + ".bedgraph"
    }
}