version development

workflow annotations{
    input{
        File vcf
        String name
        String species = "human"
        Int threads = 8
        File reference
        #Boolean offline = true
        Boolean database = false
        File ensembl_cache
        File ensembl_plugins
        String destination
    }

    call vep_annotation{
            input: vcf = vcf, ensembl_cache = ensembl_cache, name = name+"_variant_annotations.tsv", ensembl_plugins = ensembl_plugins, fasta = reference
        }

    call copy as copy_annotations{
            input:
            destination = destination + "/annotations",
            files =[
                   vep_annotation.out,
                   vep_annotation.summary
            ]
        }
}

task vep_annotation {
    input {
        File vcf
        String name = "variant_effect_output.tsv"
        String species = "homo_sapiens"
        Int threads = 8
        Boolean database = false
        File fasta
        #Boolean offline = true
        File ensembl_cache
        File ensembl_plugins
    }

    #TODO add VEP plugins http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html

    command {
        set -e
        vep --verbose --input_file ~{vcf} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything --fasta ~{fasta} \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache} --dir_plugins ~{ensembl_plugins}
    }
    #  --gene_phenotype --biotype --uniprot --symbol --allele_number --total_length --allele_number --regulatory --af

    runtime {
        docker: "ensemblorg/ensembl-vep:release_102.0"
    }

    output {
        File out = name
        File summary = name+ "_summary.html"
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