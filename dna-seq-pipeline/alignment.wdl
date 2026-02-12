version development

import "https://raw.githubusercontent.com/dna-seq/dna-seq/master/dna-seq-pipeline/clean_reads.wdl" as cleaner
import "https://raw.githubusercontent.com/dna-seq/dna-seq/master/dna-seq-pipeline/align_reads.wdl" as aligner


workflow alignment {
    input {
        Array[File]+ reads
        File reference
        File? reference_index #used only if bwa-mem2 is an aligner
        String name
        String destination
        String sequence_aligner = "bwa_mem2"
        Boolean copy_cleaned = true
        Int align_threads# = 12
        Int sort_threads# = 12
        Int max_memory_gb# = 36
        Int coverage_sampling# = 1000
        String? gencore_quality
        Boolean markdup = false
        Int compression = 9
    }

    call cleaner.clean_reads as clean_reads { input: run = name, folder = destination + "/" + "cleaned", reads = reads, copy_cleaned = copy_cleaned, is_paired = true}

    call aligner.align_reads as align_reads{
        input:
            reads = clean_reads.out.cleaned_reads,
            reference = reference,
            run = name,
            max_memory_gb = max_memory_gb,
            align_threads = align_threads,
            sort_threads = sort_threads,
            destination = destination + "/" + "aligned",
            aligner = sequence_aligner,
            markdup = markdup,
            compression = compression,
            reference_index = reference_index
    }


    output {
       CleanedRun cleaned_run =  clean_reads.out
       AlignedRun out =  align_reads.out
    }
}
#TODO: integrate tasks below somehow

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
            docker: "staphb/samtools:2.31.1" 
            maxRetries: 2
          }

    output {
        File out = name + ".bedgraph"
    }
}


task get_rg {
    input {
        Boolean rg_use_source
        File? rg_source
        String ID
        String LB
        String PL
        String PU
        String SM
    }

    command {
        touch error.txt
        if [ "~{rg_use_source}" != 'true' ]; then
        echo ~{"@RG\\\\\\\\tID:" + ID + "\\\\\\\\tLB:" + LB + "\\\\\\\\tPL:" + PL + "\\\\\\\\tPU:" + PU + "\\\\\\\\tSM:" + SM}
        else
        samtools view -H ~{rg_source} | grep '^@RG' | sed 's/'$'\t''/\\\\t/g'
        fi
    }

    runtime {
        docker_cpu: "1"
        docker: "staphb/samtools:1.23"
    }

    output {
        File e = "error.txt" #terrible hack
        String rg = read_string(stdout())
    }
}
