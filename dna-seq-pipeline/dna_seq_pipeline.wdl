version development

import "common.wdl" as common
import "quality.wdl" as quality
import "alignment.wdl" as alignment
#import "variant_calling.wdl" as variants
#import "annotations.wdl" as annotations

workflow dna_seq_pipeline {
    input {
        Array[File] reads
        String destination
        #Boolean is_paired
        File genome
        String name
    }


    call quality.quality as cleaning{
        input:
            reads = reads, is_paired = true#is_paired
    }

    call common.copy as copy_quality {
        input:
            destination = destination + "/cleaned",
            files = [cleaning.reads_cleaned[0], cleaning.reads_cleaned[1], cleaning.report_json, cleaning.report_html]
    }    

    call alignment.alignment as aligning {
        input:
            reads = cleaning.reads_cleaned,
            reference = genome,
            name = name
    }
    
    call common.copy as copy_alignemnt {
        input:
              destination = destination + "/aligned",
              files = [aligning.out, aligning.html, aligning.json]
    }

    output {
        File out = aligning.out
        File html = aligning.html
        File json = aligning.json
    }
}