version development

#Deep-variant-based pipeline
#using https://github.com/google/deepvariant/blob/r1.8/docs/deepvariant-quick-start.md as reference

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


workflow DeepVariant {
    input {
        File bam
        File bai
        String? regions
        String name = "output"
        File reference
        File reference_fai
        Int threads = 8
        String destination
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }

    call go_deep {
        input:  bam = bam, bai = bai,
            regions = regions,
            name = name,
            reference=reference, 
            reference_fai = reference_fai,
            threads = threads, 
            mode = mode
    }
    call files.copy as copy_deepvariant {
        input: files = [go_deep.vcf, go_deep.gvcf, go_deep.report, go_deep.interim], destination = destination
    }

    output {
        File vcf = copy_deepvariant.out[0]
        File gvcf = copy_deepvariant.out[1]
        File report = copy_deepvariant.out[2]
        File interim = copy_deepvariant.out[3]
    }
}
task go_deep{
    input {

        File bam
        File bai
        String name
        File reference
        File reference_fai
        Int threads
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
        String? regions
    }
    String output_vcf = name+".vcf"
    String output_gvcf = name+".g.vcf"

    command {
        ln -s ~{bam} .
        ln -s ~{bai} .
        ln -s ~{reference} .
        ln -s ~{reference_fai} .
        mkdir interim
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{mode} \
        --ref=~{basename(reference)} \
        --reads=~{basename(bam)} ~{"--regions " + regions} \
        --output_vcf=~{output_vcf} \
        --output_gvcf=~{output_gvcf} \
        --intermediate_results_dir interim \
        --num_shards=~{threads} \
        --vcf_stats_report
    }
    #    --report_title=~{name + ".visual_report.html"} \
    runtime {
        docker: "google/deepvariant:1.8.0"
    }

    output {
        File vcf = output_vcf
        File gvcf = output_gvcf
        File report = name + ".visual_report.html"
        File interim = "interim"
    }
}
