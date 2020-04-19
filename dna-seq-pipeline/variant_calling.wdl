version development

workflow variant_calling {

    input {
        File bam
        File bai
        File referenceFasta
        File referenceFai
        Int threads = 16
        Int max_memory = 28
        String name
        File ensembl_cache
        File ensembl_plugins
    }

      call smoove{
            input:
                bam = bam, bai = bai, reference = referenceFasta, reference_index = referenceFai, sample = name
        }

    call strelka2_germline{
        input:
            bam = bam, bai = bai, reference_fasta = referenceFasta, reference_fai= referenceFai, cores = threads
            #,indel_candidates = manta_germline_sv.manta_indel_candidates
    }

    call vep_annotation{
        input: vcf = strelka2_germline.variants, ensembl_cache = ensembl_cache, name = name+"_variant_annotations.tsv", ensembl_plugins = ensembl_plugins
    }

    call vep_annotation as vep_annotation_smoove{
        input: vcf = smoove.smooveVcf, ensembl_cache = ensembl_cache, name = name+"_variant_smoove_annotations.tsv", ensembl_plugins = ensembl_plugins
    }

    output {
        File variants = strelka2_germline.variants
        File variantsIndex = strelka2_germline.variantsIndex
        File variants_SV = smoove.smooveVcf
        File results_SNP = strelka2_germline.results
        File results_SV = smoove.out
        File annotations = vep_annotation.out
        File vep_summary = vep_annotation.summary
        File annotations_smoove = vep_annotation_smoove.out
    }
}


task strelka2_germline{
    input {
        String runDir = "./strelka_run"
        File bam
        File bai
        File reference_fasta
        File reference_fai
        #File? indel_candidates # ~{"--indelCandidates " + indel_candidates} \
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


task smoove {
    input {
        File bam
        File bai
        File reference
        File reference_index
        String sample
        String outputDir = "./smoove"
        Int max_memory = 16
        Int max_cores = 8
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        smoove call \
        --outdir ~{outputDir} \
        --name ~{sample} \
        --fasta ~{reference} \
        ~{bam}
    }

    output {
        File out = outputDir
        File smooveVcf = outputDir + "/" + sample + "-smoove.vcf.gz"
    }

    runtime {
        docker: "brentp/smoove@sha256:1bbf81b1c3c109e62c550783c2241acc1b10e2b161c79ee658e6abd011000c67"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{max_cores}"
    }
}




task vep_annotation {
    input {
        File vcf
        String name = "variant_effect_output.tsv"
        String species = "homo_sapiens"
        Int threads = 16
        Boolean database = false
        #Boolean offline = true
        File ensembl_cache
        File ensembl_plugins
    }

    #TODO add VEP plugins http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html

    command {
        set -e
        vep --verbose --input_file ~{vcf} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache} --dir_plugins ~{ensembl_plugins}
    }
    #  --gene_phenotype --biotype --uniprot --symbol --allele_number --total_length --allele_number --regulatory --af

    runtime {
        docker: "ensemblorg/ensembl-vep:release_99.2"
    }

    output {
        File out = name
        File summary = name+ "_summary.html"
    }
}
