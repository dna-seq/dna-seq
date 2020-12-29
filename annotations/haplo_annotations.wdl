version development

workflow haplo_annotations{
    input{
        File vcf
        String name
        String species = "homo_sapiens"
        Int threads = 16
        File reference
        #Boolean offline = true
        Boolean database = false
        File ensembl_cache
        File ensembl_plugins
        String destination
        #File disease_associations
        File G2P
        Boolean check_sorted = true
        Int buffer_size = 5000
    }

    call haplo{
        input: vcf = vcf,
            ensembl_cache = ensembl_cache,
            name = name+"_variant_annotations.tsv",
            ensembl_plugins = ensembl_plugins,
            fasta = reference,
            species = species,
        #disease_associations = disease_associations,
            G2P = G2P,
            check_sorted = check_sorted,
            buffer_size = buffer_size
    }

    call copy as copy_annotations{
        input:
            destination = destination + "/annotations",
            files =[
                   haplo.out,
                   haplo.summary
                   ]
    }
}

task haplo {
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
        #File disease_associations
        File G2P
        Boolean check_sorted
        Int buffer_size = 5000
    }

    #TODO add haplo plugins http://www.ensembl.org/info/docs/tools/haplo/script/haplo_plugins.html

    command {
        set -e
        haplo --verbose --input_file ~{vcf} -o ~{name} --tab --species ~{species} --fork ~{threads} --everything --fasta ~{fasta} \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache} --dir_plugins ~{ensembl_plugins} \
        --symbol --check_existing ~{if(check_sorted) then "" else "--no_check_variants_order"} \
        --max_sv_size 1000000000 --buffer_size ~{buffer_size}   \
        --plugin G2P,file=~{G2P},html_report=g2p_report.html,txt_report=g2p_report.txt
    }
    #--plugin DisGeNET,file=all_variant_disease_pmid_associations.tsv.gz \
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