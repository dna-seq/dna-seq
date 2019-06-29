version development

import "common.wdl" as common
import "quality.wdl" as quality
import "alignment.wdl" as alignment
import "variant_calling.wdl" as variants
import "annotations.wdl" as annotations

workflow dna_seq_pipeline {
    input {
        Array[File] reads
        Boolean is_paired
        File genome
        File gtf
        String results
    }
}