version development

workflow alignment {
    input {
        Array[File] reads
        File reference
        String name
        Int threads = 4
        String? gencore_quality
    }

    call minimap2 {
        input:
            reads = reads,
            reference = reference,
            name = name,
            threads = threads
    }

    call samtools_conversion {
        input:
            sam = minimap2.out
    }

    call samtools_sort {
        input:
            bam = samtools_conversion.out,
            threads = 4
    }

    call gencore{
        input:
            reference = reference,
            sorted_bam = samtools_sort.out,
            name = name,
            quality  = gencore_quality
    }

    call samtools_sort as sort_gencore{
        input:
            bam = gencore.out,
            threads = 4
    }

    output {
       File bam = sort_gencore.out
       File bai = sort_gencore.bai
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
    }

    command {
        minimap2 -ax sr  -t ~{threads} -2 ~{reference} ~{sep=' ' reads} > ~{name}.sam
    }

    runtime {
        docker: "quay.io/biocontainers/minimap2@sha256:e7ec93ae6c9dcdad362333932ad8f509cbee6de691c72cb9e600972c770340f4" #2.17--h8b12597_1
        maxRetries: 2
      }

    output {
      File out = name + ".sam"
    }
}

task samtools_conversion {
    input {
        File sam
    }

    String name = basename(sam, ".sam")

    command {
       samtools view -bS ~{sam} > ~{name}.bam
    }

    runtime {
        docker: "quay.io/biocontainers/samtools@sha256:70581cfc34eb40cb9b55e49cf5805fce820ec059d7bca9bbb762368ac3c1ac0a" #1.10--h9402c20_2
        maxRetries: 2
      }

    output {
        File out = name + ".bam"
      }
}


task samtools_sort {
    input {
        File bam
        Int threads
    }

    String name = basename(bam, ".bam")

    command {
       samtools sort ~{bam} --threads ~{threads} -o ~{name}_sorted.bam
       samtools index ~{name}_sorted.bam  ~{name}_sorted.bai
    }

    runtime {
        docker: "biocontainers/samtools@sha256:da61624fda230e94867c9429ca1112e1e77c24e500b52dfc84eaf2f5820b4a2a" #v1.9-4-deb_cv1
        maxRetries: 2
      }

    output {
        File out = name + "_sorted.bam"
        File bai = name + "_sorted.bai"
      }
}

task gencore {
    input {
        File reference
        File sorted_bam
        String name
        Int supporting_reads = 1
        Float ratio_threshold = 0.8
        String? quality #"--high_qual"
    }
    command {
        gencore --ratio_threshold=~{ratio_threshold} -s ~{supporting_reads} ~{quality} -i ~{sorted_bam} -o ~{name}.bam -r ~{reference}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/gencore@sha256:2bde467b5d57ee17397a9a2aeb82e18279032f12d760baaad9e6d38227157565"
    }

    output {
        File out = name + ".bam"
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