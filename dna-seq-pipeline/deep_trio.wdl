version development

#alternative deep-variant-based pipeline, work in progress
#using https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-quick-start.md as reference
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow DeepTrio {
    input {
        File parent1_bam
        File parent1_bai
        File parent2_bam
        File parent2_bai
        File child_bam
        File child_bai
        String parent1_sample_name
        String parent2_sample_name
        String child_sample_name
        File reference
        File reference_fai
        Int threads = 8
        String destination
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }

    call go_deep {
        input:
            parent1_bam = parent1_bam, parent1_bai = parent1_bai,
            parent2_bam = parent2_bam, parent2_bai = parent2_bai,
            child_bam = child_bam, child_bai = child_bai,
            parent1_name = parent1_sample_name,
            parent2_name = parent2_sample_name,
            child_name = child_sample_name,
            reference=reference, 
            reference_fai = reference_fai,
            threads = threads, 
            mode = mode
    }
    call files.copy as copy_deepvariant {
        input:
            files = [ go_deep.child_vcf, go_deep.child_gvcf, go_deep.child_report,
                    go_deep.parent1_vcf, go_deep.parent1_gvcf, go_deep.parent1_report,
                    go_deep.parent2_vcf, go_deep.parent2_gvcf, go_deep.parent2_report
                    ], destination = destination
    }

    output {
        File child_vcf = copy_deepvariant.out[0]
        File child_gvcf = copy_deepvariant.out[1]
        File child_report = copy_deepvariant.out[2]
		File parent1_vcf = copy_deepvariant.out[3]
        File parent1_gvcf = copy_deepvariant.out[4]
        File parent1_report = copy_deepvariant.out[5]  
		File parent2_vcf = copy_deepvariant.out[6]
        File parent2_gvcf = copy_deepvariant.out[7]     
		File parent2_report = copy_deepvariant.out[8]  
    }
}
task go_deep{
    input {
        File parent1_bam
        File parent1_bai
        File parent2_bam
        File parent2_bai
        File child_bam
        File child_bai
        String parent1_name
        String parent2_name
        String child_name
        File reference
        File reference_fai
        Int threads
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }
    String child_output_vcf = child_name+".vcf"
    String child_output_gvcf = child_name+".g.vcf"
	String child_output_visual = child_name+".visual_report.html"
    String parent1_output_vcf = parent1_name+".vcf"
    String parent1_output_gvcf = parent1_name+".g.vcf"
	String parent1_output_visual = parent1_name+".visual_report.html"
    String parent2_output_vcf = parent2_name+".vcf"
    String parent2_output_gvcf = parent2_name+".g.vcf"
	String parent2_output_visual = parent2_name+".visual_report.html"

    command {
        ln -s ~{parent1_bam} .
        ln -s ~{parent1_bai} .
        ln -s ~{parent2_bam} .
        ln -s ~{parent2_bai} .
        ln -s ~{child_bam} .
        ln -s ~{child_bai} .
        ln -s ~{reference} .
        ln -s ~{reference_fai} .
        mkdir -p logs
        mkdir -p interim
        
        /opt/deepvariant/bin/deeptrio/run_deeptrio \
        --model_type=~{mode} \
        --num_shards=~{threads} \
        --logging_dir logs \
        --intermediate_results_dir interim \
        --ref=~{basename(reference)} \
        --reads_child=~{basename(child_bam)} \
        --reads_parent1=~{basename(parent1_bam)} \
        --reads_parent2=~{basename(parent2_bam)} \
        --sample_name_child=~{child_name} \
        --sample_name_parent1=~{parent1_name} \
        --sample_name_parent2=~{parent2_name} \
        --output_vcf_child ~{child_output_vcf} \
        --output_vcf_parent1 ~{parent1_output_vcf} \
        --output_vcf_parent2 ~{parent2_output_vcf} \
        --output_gvcf_child ~{child_output_gvcf} \
        --output_gvcf_parent1 ~{parent1_output_gvcf} \
        --output_gvcf_parent2 ~{parent2_output_gvcf} \
        --runtime_report \
        --vcf_stats_report

    }

    runtime {
        docker: "google/deepvariant:deeptrio-1.5.0"
    }

    output {
        File child_vcf = child_output_vcf
        File child_gvcf = child_output_gvcf
        File child_report = child_output_visual
		File parent1_vcf = parent1_output_vcf
        File parent1_gvcf = parent1_output_gvcf
        File parent1_report = parent1_output_visual
		File parent2_vcf = parent2_output_vcf
        File parent2_gvcf = parent2_output_gvcf
        File parent2_report = parent2_output_visual
    }
}