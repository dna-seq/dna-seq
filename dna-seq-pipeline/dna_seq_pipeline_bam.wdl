version 1.0

import "common.wdl" as common
import "quality.wdl" as quality
import "alignment.wdl" as alignment
#import "variant_calling.wdl" as variants
#import "annotations.wdl" as annotations

workflow dna_seq_pipeline_bam {

    input {
        File bam
        String destination
        #Boolean is_paired
        File genome
        String name
    }

    call alignment.alignment {

    }

    output {
        File out = alignment.out
        File html = aligning.html
        File json = aligning.json
    }
}