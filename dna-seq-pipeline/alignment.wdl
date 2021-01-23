version development

workflow alignment {
    input {
        Array[File]+ reads
        File reference
        String name
        String destination
        Int align_threads# = 12
        Int sort_threads# = 12
        Int max_memory_gb# = 36
        Int coverage_sampling# = 1000
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
            quality = gencore_quality,
            coverage_sampling = coverage_sampling,
            max_memory = max_memory_gb
    }


    call copy as copy_alignment {
        input:
            destination = destination,
            files = [gencore.bam, gencore.bai, gencore.html, gencore.json]
    }


    output {
       File bam =  copy_alignment.out[0]
       File bai = copy_alignment.out[1]
       File html = copy_alignment.out[2]
       File json = copy_alignment.out[3]
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
        docker: "quay.io/comp-bio-aging/minimap2@sha256:69e9515a0cb5b5e9f47c3d0f95700d5064f0db04f86d49b6626b66a012daf0a5" #latest
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
    }

    runtime {
        docker: "quay.io/biocontainers/sambamba@sha256:ae92faef4c53a632b2120dfffa7b6dcfe5366a0647e61bbbd6188aedc89da4e8" #:0.8.0--h984e79f_0
        maxRetries: 1
        docker_memory: "~{gb_per_thread * (threads+1)}G"
        docker_cpu: "~{threads+1}"
        docker_swap: "~{gb_per_thread * (threads+1) * 2}G"
      }

    output {
      File out = name + ".sorted.bam"
      File bai = name + ".sorted.bam.bai"
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
        Int coverage_sampling# = 1000
    }
    command {
        gencore --coverage_sampling ~{coverage_sampling} --ratio_threshold=~{ratio_threshold} -s ~{supporting_reads} ~{quality} -i ~{sorted_bam} -o ~{name}.bam -r ~{reference}
        samtools index ~{name}.bam  ~{name}.bam.bai
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker: "quay.io/comp-bio-aging/gencore@sha256:14b0da6c870766e04ea80a3d010ee593bf6a0bd071c5d4cdee002095a632a828"
    }

    output {
        File bam = name + ".bam"
        File bai = name + ".bam.bai"
        File html = "gencore.html"
        File json = "gencore.json"
    }
}

task coverage {
    input {
        File bam
        #takes a lot of time and consumes space, so far it is switched off
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
