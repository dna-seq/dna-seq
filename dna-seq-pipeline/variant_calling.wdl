version development

workflow variant_calling {

    input {
        File bam
        File bai
        File referenceFasta
        File referenceFai
        Int threads = 16
        Int max_memory = 28
    }

    call strelka2_germline{
        input:
            bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai, cores = threads
    }

    call parlament2{
        input:
            bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai
    }

    call vep_annotation{
        input: vcf = strelka2_germline.variants
    }

    output {
        File variants = strelka2_germline.variants
        File variantsIndex = strelka2_germline.variantsIndex
        File results_SNP = strelka2_germline.results
        File results_CNV = parlament2.out
        File annotations = vep_annotation.out
        File vep_summary = vep_annotation.summary
    }
}

task strelka2_germline{
    input {
        String runDir = "./strelka_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        Boolean exome = false
        Boolean rna = false

        Int cores = 8
        Int max_memory = 28
    }

    command {
        configureStrelkaGermlineWorkflow.py \
        --bam ~{bam} \
        --ref ~{reference_fasta} \
        --runDir ~{runDir} \
        ~{true="--rna" false="" rna}
        ~{true="--exome" false="" exome} \
        ~{runDir}/runWorkflow.py \
        -m local \
        -j ~{cores} \
        -g ~{max_memory}
    }

    output {
        File variants = runDir + "/results/variants/variants.vcf.gz"
        File variantsIndex = runDir + "/results/variants/variants.vcf.gz.tbi"
        File results = runDir
    }

    runtime {
        docker: "quay.io/biocontainers/strelka:2.9.10--0"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{cores}"
    }
}

# This task takes a BAM file and runs parliament2 on it.
task parlament2 {

  input
  {
    File bam
    File bai
    File reference_fai
    File reference_fasta
  }

  String name = basename(bam, ".bam")

  command {
      /home/dnanexus/parliament2.sh --bam ~{bam} --bai ~{bai} --fai ~{reference_fai} -r ~{reference_fasta} \
      --breakdancer --breakseq --manta --cnvnator --lumpy \
      --delly_deletion --delly_insertion --delly_inversion --delly_duplication \
      --genotype --svviz
  }

  runtime {
    docker: "dnanexus/parliament2:latest"
  }

  output {
    File out = "output"
    #Array[File] vcfs = glob("output/*.vcf")
    #Array[File] sv_caller_results = glob("output/sv_caller_results/*")
    #Array[File] svtyped_vcfs = glob("output/svtyped_vcfs/*.vcf")
  }

}


task vep_annotation {
    input {
        File vcf
        String name = "variant_effect_output.tsv"
        String species = "human"
        Int threads = 8
        Boolean offline = true
    }

    #TODO add VEP plugins http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html

    command {
        vep \
        --input_file ~{vcf} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything ~{if(offline) then "--offline" else ""} \
        --gene_phenotype --biotype --uniprot --symbol --allele_number --total_length --allele_number --regulatory --af
    }

    runtime {
        docker: "ensemblorg/ensembl-vep:release_99.2"
    }

    output {
        File out = name
        File summary = "variant_effect_output.txt_summary.html"
    }
}