version development

workflow alignment {
    input {
        Array[File] reads
        File reference
        String name
        Int threads = 4
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
            bam = samtools_conversion.out
    }

    call gencore{
        input:
            reference = reference,
            sorted_bam = samtools_sort.out,
            name = name
    }

    output {
       File out = gencore.out
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
        docker: "quay.io/biocontainers/samtools@sha256:f8d74e566c4a250d35bc3194298a034ae096e348532188ca743354a7d5fff50c" #1.9--h10a08f8_12
        maxRetries: 2
      }

    output {
        File out = name + ".bam"
      }
}


task samtools_sort {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
       samtools sort ~{bam}  -o ~{name}_sorted.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
        maxRetries: 2
      }

    output {
        File out = name + "_sorted.bam"
      }
}

task gencore {
    input {
        File reference
        File sorted_bam
        String name
        Int supporting_reads = 2
        Float ratio_threshold = 0.8
        String? quality #"--high_qual"
    }
    command {
        gencore --ratio_threshold=~{ratio_threshold} -s ~{supporting_reads} ~{quality} -i ~{sorted_bam} -o ~{name}.bam -r ~{reference} --coverage_sampling=50000
    }

    runtime {
        docker: "quay.io/comp-bio-aging/gencore"
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
            docker: "quay.io/biocontainers/bedtools@sha256:a0bb135afdec53be4b953a9a8efbc801cdb90706e6e63e11e3f60b06b8444f78" #2.23.0--he941832_1
            maxRetries: 2
          }

    output {
        File out = name + ".bedgraph"
    }
}