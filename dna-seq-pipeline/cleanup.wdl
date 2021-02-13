version development

workflow cleanup {
    input {
        Array[File]+ reads
        String destination
        Int threads# = 16
        Boolean is_paired = true #now it supports only paired, TODO: fix for single
    }

    call fastp { input: reads = reads, threads = threads, is_paired = is_paired }
    call copy as copy_cleaned { input: destination = destination, files = fastp.out }


    output {
        Array[File]+ reads_cleaned = [copy_cleaned.out[0], copy_cleaned.out[1]]
        File report_json = copy_cleaned.out[2]
        File report_html = copy_cleaned.out[3]
    }

}


task fastp {
    input {
        Array[File]+ reads
        Int threads
        Boolean is_paired
    }

    command {
        fastp --cut_front --cut_tail --cut_right --overrepresentation_analysis \
        -i ~{reads[0]} -w ~{threads} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
        ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
                                    then [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz", basename(reads[1], ".fastq.gz") + "_cleaned.fastq.gz"]
                                    else [basename(reads[0], ".fastq.gz") + "_cleaned.fastq.gz"]
        Array[File] out = if( is_paired )
                          then [reads_cleaned[0], reads_cleaned[1], report_json, report_html] else
                          [reads_cleaned[0], report_json, report_html]
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
