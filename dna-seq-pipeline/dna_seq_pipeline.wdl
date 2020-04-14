version development

import "alignment.wdl" as alignment
import "recalibration.wdl" as recalibration
import "variant_calling.wdl" as variant_calling

struct Reference{
    File genome #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa
    File fai #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
    File dict #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.dict
}

workflow dna_seq_pipeline {
    input {
        Array[File]+ reads
        String destination
        Boolean is_paired = true
        Reference reference
        String name
        Int threads = 8
    }

    call fastp { input: reads = reads, is_paired = is_paired }
    call copy as copy_cleaned { input: destination = destination + "/fastq/cleaned", files = fastp.out }

    call alignment.alignment as align{
        input:
          reads = fastp.reads_cleaned,
          reference = reference.genome,
          name = "",
          threads = threads
    }

    call copy as copy_alignment{
        input:
        destination = destination + "/bam/aligned",
        files = [ align.bam, align.bai, align.html, align.json]
    }

    call recalibration.recalibration as recalib {
        input:
            reference = reference.genome,
            referenceDict = reference.dict,
            referenceFai = reference.fai,
            bam = copy_alignment.out[0],
            bai = copy_alignment.out[1]
    }

    call copy as copy_recalibration{
        input:
        destination = destination + "/bam/recalibrated",
        files = [ recalib.recalibratedBam, recalib.recalibratedBamIndex, recalib.recalibratedBamMd5, recalib.report]
    }

    output {
        File recalibratedBam = copy_recalibration.out[0]
        File recalibratedBamIndex =  copy_recalibration.out[1]
        File recalibratedBamMd5 =  copy_recalibration.out[2]
        File report = copy_recalibration.out[3]
        File results_folder = destination
    }


}

task fastp {
    input {
        Array[File]+ reads
        Boolean is_paired
    }

    command {
        fastp --cut_front --cut_tail --cut_right --overrepresentation_analysis \
            -i ~{reads[0]} -o ~{basename(reads[0], ".fastq.gz")}_cleaned.fastq.gz \
            ~{if( is_paired ) then "--detect_adapter_for_pe " + "--correction -I "+reads[1]+" -O " + basename(reads[1], ".fastq.gz") +"_cleaned.fastq.gz" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:ac9027b8a8667e80cc1661899fb7e233143b6d1727d783541d6e0efffbb9594e" #0.20.0--hdbcaa40_0
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