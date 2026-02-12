version development

# production version
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#local debug version (uncomment for debugging and comment the production version)
#import "../common/files.wdl" as files

struct AlignedRun {
    String run
    String folder
    File bam
    File bai
    File flagstat
    String aligner
}

workflow align_reads {

    input {
        Array[File] reads
        File reference
        File? reference_index #if no index is provided that computes it inplace
        String run
        Int max_memory_gb = 42
        Int align_threads = 12
        Int sort_threads = 12
        Int gb_per_thread = 3
        Boolean markdup = false
        String destination
        String aligner = "bwa-mem2"
        Boolean markdup = false
        Int compression = 9
        Boolean readgroup = true
    }

    Boolean is_bwa = (aligner == "bwa_mem2" || aligner == "bwa-mem2")

    if(aligner == "minimap2") {
        call minimap2 {
            input:
                reads = reads,
                reference = reference,
                name = run,
                threads = align_threads,
                max_memory = max_memory_gb,
                readgroup = readgroup,
                compression = compression
        }
    }
    if(is_bwa) {
           call bwa_mem2 {
               input:
                    reads = reads,
                    reference = reference,
                    reference_index = reference_index,
                    name = run,
                    threads = align_threads,
                    max_memory = max_memory_gb,
                    readgroup = readgroup,
                    compression = compression
           }
     }

    File unsorted_bam = select_first([if(is_bwa) then bwa_mem2.bam else minimap2.bam])

    call sambamba_sort {
        input:
            unsorted_bam = unsorted_bam,
            threads = sort_threads,
            gb_per_thread = gb_per_thread,
            markdup = markdup,
            compression = compression
    }

    call files.copy as copy_sorted_bam{
        input:
            destination = destination,
            files = [sambamba_sort.sorted_bam, sambamba_sort.sorted_bai, sambamba_sort.flagstat]
    }

    output {
       AlignedRun out = object {run: run, folder: destination, bam: copy_sorted_bam.out[0], bai: copy_sorted_bam.out[1], flagstat: copy_sorted_bam.out[2], aligner: aligner}
    }
}

task minimap2 {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
        Int max_memory
        Int compression = 9
        Boolean readgroup = true
    }

    #TODO for readgroups investigate https://angus.readthedocs.io/en/2017/Read_group_info.html

    command {
        minimap2 ~{if(readgroup) then "-R '@RG\tID:" +name +"'" else ""} -ax map-ont -t ~{threads} -2 ~{reference} ~{sep=' ' reads} \
        | sambamba view -t ~{threads} -l ~{compression} -S -f bam -o ~{name}.bam /dev/stdin
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker: "justdnaseq/minimap2:2.30"
        maxRetries: 2
    }

    output {
        File bam = name + ".bam"
    }
}


task bwa_mem2 {
    input {
        Array[File] reads
        File reference
        File? reference_index
        String name
        Int threads
        Int max_memory
        Int compression = 9
        Boolean readgroup = true
    }

    String ref_name = basename(reference)
    Boolean has_index = defined(reference_index)

    command {
        ln -s ~{reference} ~{ref_name}
        ~{if(has_index) then "ln -s " + reference_index + " " + basename(select_first([reference_index])) else "bwa-mem2 index " +  ref_name}
        bwa-mem2 mem ~{if(readgroup) then "-R '@RG\tID:" +name +"'" else ""} -t ~{threads} ~{if(has_index) then basename(select_first([reference_index]))+"/"+ ref_name else ref_name} ~{sep=' ' reads} \
        | sambamba view -t ~{threads} -l ~{compression} -S -f bam -o ~{name}.bam /dev/stdin
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker: "justdnaseq/bwa-mem2:2.3"
        maxRetries: 2
    }

    output {
        File bam = name + ".bam"
    }
}


task sambamba_sort {
    input {
        File unsorted_bam
        Int threads
        Int gb_per_thread = 3
        Int compression = 9
        Boolean markdup = false
    }

    String name = basename(unsorted_bam, ".bam")
    String suffix = if(markdup) then ".sorted_markdup.bam" else ".sorted.bam"

    command {
        ln -s ~{unsorted_bam} ~{basename(unsorted_bam)}
        sambamba sort -l ~{compression} -m ~{gb_per_thread}G -t ~{threads} -p ~{basename(unsorted_bam)}
        ~{if(markdup) then "sambamba markdup -t " + threads + " -l " + compression +" -p " + name + ".sorted.bam " + name + suffix else ""}
        sambamba flagstat -t ~{threads} -p ~{name + suffix} > ~{name + suffix + ".flagstat"}
    }

    runtime {
        docker: "justdnaseq/sambamba:1.0.1"
        maxRetries: 1
        docker_memory: "~{gb_per_thread * (threads+1)}G"
        docker_cpu: "~{threads+1}"
        docker_swap: "~{gb_per_thread * (threads+1) * 2}G"
    }

    output {
        File sorted_bam = name + suffix
        File sorted_bai = name + suffix + ".bai"
        File flagstat = name + suffix + ".flagstat"
    }
}
