version development


workflow dna_seq_pipeline_fastq {
    input {
        Array[File]+ reads
        String destination
        #Boolean is_paired
        File genome
        String name
        Int threads = 3
    }

    call fastp as cleaning { input: reads = reads,  adapters = ["AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"], q = 32 }

    call copy as copy_quality {
        input:
            destination = destination + "/cleaned",
            files = [cleaning.reads_cleaned[0], cleaning.reads_cleaned[1], cleaning.report_json, cleaning.report_html]
    }

    call minimap2 {
        input:
            reads = cleaning.reads_cleaned,
            reference = genome,
            name = name,
            threads = 4
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
            reference = genome,
            sorted_bam = samtools_sort.out,
            name = name
    }

    call samtools_sort as sort_gencore{
        input:
            bam = gencore.out,
            threads = 4
    }
    call copy as copy_alignemnt {
        input:
              destination = destination + "/aligned",
              files = [sort_gencore.out, gencore.html, gencore.json]
    }


    output {
        File out =  sort_gencore.out
        File html = gencore.html
        File json = gencore.json
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
        docker: "quay.io/comp-bio-aging/gencore@sha256:71d64d00e1a50478136ed292d2ca28a1b4ae85856b61d1282989b47ed22e5d2e"
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


task fastp {
    input {
        Array[File] reads
        Array[String] adapters
        Int q = 35
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis --correction \
            --adapter_sequence ~{adapters[0]} --adapter_sequence_r2 ~{adapters[1]} -q ~{q} \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fq.gz")}_cleaned.fq.gz \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq.gz") +"_cleaned.fq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
            then [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz", basename(reads[1], ".fq.gz") + "_cleaned.fq.gz"]
            else [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz"]
    }
}


task quast {

    input {
        Array[File]+ contigs
        File? reference
        Int? threads = 4
        File? features
        String output_folder = "results"
        Int min_contig = 50
        String? type
    }

    command {
        quast.py ~{if defined(reference) then "--reference " + reference else ""} \
         ~{if defined(threads) then "--threads " + threads else ""} ~{sep=" " contigs} \
         --output ~{output_folder} \
         ~{if defined(features) then "--features " + features + (if(defined(type)) then "--type " + type else "") else "" } \
         --min-contig ~{min_contig} \
         ~{sep=" " contigs}
    }

    runtime {
        docker: "quay.io/biocontainers/quast@sha256:89c337541c3bc92bed901b6215231a5b6f18bed86e25b5f94a97fee73d0e7c13" #5.0.2--py27pl526ha92aebf_0
    }

    output {
        File out = output_folder
    }
}
