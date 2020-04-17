version development

import "alignment.wdl" as alignment
import "variant_calling.wdl" as vc

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
        Int align_threads = 8
        Int sort_threads = 12
        Int variant_calling_threads = 8
        String name
    }

    call fastp { input: reads = reads, is_paired = is_paired }
    call copy as copy_cleaned { input: destination = destination + "/fastq/cleaned", files = fastp.out }

    call alignment.alignment as align{
        input:
          reads = fastp.reads_cleaned,
          reference = reference.genome,
          name = name,
          align_threads = align_threads,
          sort_threads = sort_threads
    }

    call copy as copy_alignment{
        input:
        destination = destination + "/bam/aligned",
        files = [ align.bam, align.bai, align.html, align.json]
    }

    call vc.variant_calling as variant_calling{
        input:
                bam = copy_alignment.out[0],
                bai = copy_alignment.out[1],
                referenceFasta = reference.genome,
                referenceFai = reference.fai,
                threads = variant_calling_threads
    }

    call copy as copy_variants{
        input:
        destination = destination + "/variants",
        files =[
                variant_calling.results_SNP
        ]
    }


    call copy as copy_annotations{
        input:
        destination = destination + "/variants/annotations",
        files =[
               variant_calling.annotations,
               variant_calling.vep_summary
        ]
    }

    call copy as copy_CNV{
        input:
        destination = destination + "/CNV",
        files =[
                variant_calling.results_CNV,
        ]
    }

    output {
        File results_SNP = copy_variants.out[0]
        File results_CNV =  copy_variants.out[0]
        File annotations = copy_annotations.out[0]
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