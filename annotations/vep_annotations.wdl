version development

workflow vep_annotations{
    input{
        File vcf
        File? vcf_tbi
        String name
        String species = "homo_sapiens"
        Int threads
        File reference
        #Boolean offline = true
        Boolean database = false
        File ensembl_cache
        File ensembl_plugins
        String destination
        Boolean check_sorted = true
        Int buffer_size = 5000
        File G2P
        File disease_associations
        File? disease_associations_tbi
        File clinvar
        File? clinvar_tbi
        Boolean condel = false
        Boolean conservation = false
    }

    call vep{
            input: vcf = vcf,
                ensembl_cache = ensembl_cache,
                name = name+"_variant_annotations.tsv",
                ensembl_plugins = ensembl_plugins,
                species = species,
                threads = threads,
                database = database,
                fasta = reference,
                species = species,
                disease_associations = disease_associations,
                disease_associations_tbi = disease_associations_tbi,
                G2P = G2P,
                check_sorted = check_sorted,
                buffer_size = buffer_size,
                clinvar = clinvar,
                clinvar_tbi = clinvar_tbi,
                conservation = conservation,
                condel = condel
        }

    call copy as copy_annotations{
            input:
            destination = destination + "/annotations",
            files =[
                   vep.out,
                   vep.summary,
                   vep.g2p_report
            ]
        }
}

task vep {
    input {
        File vcf
        File? vcf_tbi
        String name = "variant_effect_output.tsv"
        String species 
        Int threads
        Boolean database
        File fasta
        #Boolean offline = true
        File ensembl_cache
        File ensembl_plugins
        Boolean check_sorted
        Int buffer_size = 5000
        File G2P
        File disease_associations
        File? disease_associations_tbi
        File clinvar
        File? clinvar_tbi
        Boolean conservation
        Boolean condel
    }

    #TODO add VEP plugins http://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html

    command {
        set -e
        ln -s ~{vcf} .
        ~{"ln -s "+ vcf_tbi  +" ."}
        ~{"ln -s "+ clinvar  +" ."}
        ~{"ln -s "+ clinvar_tbi  +" ."}
        ~{"ln -s "+ disease_associations  +" ."}
        ~{"ln -s "+ disease_associations_tbi  +" ."}
        vep --verbose --input_file ~{basename(vcf)} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything --fasta ~{fasta} \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache} --dir_plugins ~{ensembl_plugins} \
        --symbol --check_existing ~{if(check_sorted) then "" else "--no_check_variants_order"} \
        --max_sv_size 1000000000 --buffer_size ~{buffer_size}   \
        ~{"--plugin G2P,file="+G2P+",html_report=g2p_report.html,txt_report=g2p_report.txt"} \
        ~{"--plugin DisGeNET,file=" + basename(disease_associations) + ",disease=1"} \
        ~{"--custom " + basename(clinvar) + ",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN"} \
        ~{if(condel) then "--plugin Condel --plugin ExACpLI" else ""} \
        ~{if(conservation) then "--plugin Conservation,method_link_type=GERP_CONSERVATION_SCORE,species_set=mammals" else ""}
    }

    #--plugin Phenotypes,output_format=json,phenotype_feature=1 \

    runtime {
        docker: "quay.io/comp-bio-aging/vep" #based on ensemblorg/ensembl-vep:release_102.0
    }

    output {
        File out = name
        File summary = name+ "_summary.html"
        File g2p_report = "g2p_report.html"
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
