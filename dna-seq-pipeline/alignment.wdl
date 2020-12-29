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
        docker: "quay.io/comp-bio-aging/minimap2@sha256:7d7ee2187e62720c88ae8facf9a0db87a5cf0517749dead432e92ec81bd30bb1" #latest
        #docker: "quay.io/biocontainers/minimap2@sha256:7f95eecc8eeee8ef8ae7e24d1d1a49ee056fb21d72aea4e2def97484f8a206c5" #2.17--hed695b0_3
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
        docker: "quay.io/biocontainers/sambamba@sha256:9ec72d3d0991c4209830e4ff17937986808c64c430780071559e7072e8317ab3" #:0.7.1--h984e79f_3
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
        docker: "quay.io/comp-bio-aging/gencore@sha256:14b0da6c870766e04ea80a3d010ee593bf6a0bd071c5d4cdee002095a632a828"
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