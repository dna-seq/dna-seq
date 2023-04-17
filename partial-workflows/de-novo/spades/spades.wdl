version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow plasmid_spades {

    input {
        Array[File] reads
        String destination
    }

    call fastp { input: reads = reads }

    call files.copy as copy_report {
        input:
            destination = destination + "/report",
            files = [fastp.report_json, fastp.report_html]
    }

    call spades {
        input: reads = reads,
    }

    call files.copy as copy_results {
        input:
            destination = destination,
            files = [spades.out]
    }

    output {
        File out = spades.out
    }
}

task fastp {
    input {
        Array[File] reads
    }

    Boolean is_paired = if(length(reads) > 1) then true else false

    command {
        fastp --cut_front --cut_tail --cut_right --trim_poly_g --trim_poly_x --overrepresentation_analysis \
        -i ~{reads[0]} -o ~{basename(reads[0], ".fq.gz")}_cleaned.fq.gz \
        ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fq.gz") +"_cleaned.fq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5" #0.20.1--h8b12597_0
    }

    output {
        File report_json = "fastp.json"
        File report_html = "fastp.html"
        Array[File] reads_cleaned = if( is_paired )
                                    then [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz", basename(reads[1], ".fq.gz") + "_cleaned.fq.gz"]
                                    else [basename(reads[0], ".fq.gz") + "_cleaned.fq.gz"]
    }
}

task spades {
    input {
        String results = "results"
        Array[File] reads
        String cut_off = "auto"
    }
    command {
        plasmidspades.py -1 ~{reads[0]} -2 ~{reads[1]} --cov-cutoff ~{cut_off} -o ~{results}
    }

    runtime {
        docker: "quay.io/biocontainers/spades@sha256:ce6565d2e86ca1a0baacca9888837c246ec1cab0f2821ed3ba40acc6e0dd5d0e" #:3.15.5--h95f258a_1"
    }

    output {
        File out = "results"
    }
}