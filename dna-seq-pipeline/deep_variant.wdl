version development

#alternative deep-variant-based pipeline, work in progress
#using https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-quick-start.md as reference

workflow DeepVariant {
    input {
        File bam
        File bai
        String? regions
        String name = "output"
        File reference #/input/ucsc.hg19.chr20.unittest.fasta
        File reference_fai
        Int threads
        String destination
        String mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }

    call go_deep {
        input:  bam = bam, bai = bai,
            regions = regions,
            name = name,
            reference=reference, reference_fai = reference_fai,
            threads = threads, mode = mode
    }
    call copy {
        input: files = [go_deep.vcf, go_deep.gvcf, go_deep.report, go_deep.interim], destination = destination
    }

    output {
        File vcf = copy.out[0]
        File gvcf = copy.out[1]
        File report = copy.out[2]
        File interim = copy.out[3]
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
        --num_shards=~{threads}
    }

    runtime {
        docker: "google/deepvariant:1.1.0"
    }

    output {
        File vcf = output_vcf
        File gvcf = output_gvcf
        File report = name + ".visual_report.html"
        File interim = "interim"
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