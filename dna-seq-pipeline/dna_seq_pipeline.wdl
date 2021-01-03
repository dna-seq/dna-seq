version development

import "cleanup.wdl" as clean
import "alignment.wdl" as alignment
import "simple_variant_calling.wdl" as vc
import "deep_variant.wdl" as deep

workflow dna_seq_pipeline {
    input {
        Array[File]+ reads
        String destination
        Boolean is_paired = true
        File reference #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa
        File reference_fai #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
        Int align_threads = 12
        Int sort_threads = 12
        Int variant_calling_threads = 16
        Int coverage_sampling = 1000
        Int max_memory_gb = 42
        String name
    }

    call clean.cleanup as cleanup{
        input:
            reads = reads, destination = destination + "/fastq/cleaned"
    }


    call alignment.alignment as align{
        input:
          reads = cleanup.reads_cleaned,
          reference = reference,
          destination = destination + "/aligned",
          name = name,
          align_threads = align_threads,
          sort_threads = sort_threads,
          coverage_sampling = coverage_sampling,
          max_memory_gb = max_memory_gb
    }


    call vc.simple_variant_calling as simple_variant_calling{
        input:
                bam = align.bam,
                bai = align.bai,
                destination = destination + "/variants",
                referenceFasta = reference,
                referenceFai = reference_fai,
                threads = variant_calling_threads,
                name = name
    }

    call deep.DeepVariant as deepvariant {
        input:
            bam = align.bam,
            bai = align.bai,
            name = name,
            reference = reference,
            reference_fai = reference_fai,
            threads = variant_calling_threads,
            destination = destination + "/variants/deepvariant",
            mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }



    output {
        File alignment = align.bam
        File deep_SNP = deepvariant.vcf
        File deep_g_SNP = deepvariant.gvcf
        File results_SNP = simple_variant_calling.results_SNP
        File results_SV =  simple_variant_calling.variants_SV
    }


}