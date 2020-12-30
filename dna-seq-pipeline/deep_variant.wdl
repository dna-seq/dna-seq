version development

#alternative deep-variant-based pipeline, work in progress
#using https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-quick-start.md as reference

workflow DeepVariant {
    input {
        File bam = "NA12878_S1.chr20.10_10p1mb.bam"
        String regions = "chr20:10,000,000-10,010,000"
        String output_vcf = "output.vcf.gz"
        String output_gvcf= "output.g.vcf.gz"
        String intermediate_results_dir
        File reference #/input/ucsc.hg19.chr20.unittest.fasta
        Int threads
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }

    call go_deep {
        input:  bam = bam, regions = regions, output_vcf=output_vcf, output_gvcf=output_gvcf, reference=reference, threads = threads, mode = mode
    }

    output {
        File vcf = go_deep.vcf
        File gvcf = go_deep.gvcf
    }
}
task go_deep{
    input {

        File bam
        String regions
        String output_vcf
        String output_gvcf
        File reference
        Int threads
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }
    command {
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{mode} \
        --reads=~{bam} \
        --regions ~{regions} \
        --output_vcf=~{output_vcf} \
        --output_gvcf=~{output_gvcf} \
        --num_shards=~{threads}
    }

    runtime {
        docker: "google/deepvariant"
    }

    output {
        File vcf = output_vcf
        File gvcf = output_gvcf
    }
}