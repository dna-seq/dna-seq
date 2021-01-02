version development

workflow haplo_annotations{
    input{
        File vcf
        File? vcf_tbi
        String name
        String species = "homo_sapiens"
        Int threads = 16
        File reference
        #Boolean offline = true
        Boolean database = false
        File ensembl_cache
        String destination
    }

    call haplo{
        input: vcf = vcf,vcf_tbi = vcf_tbi,
            ensembl_cache = ensembl_cache,
            name = name+"_variant_annotations.tsv",
            fasta = reference,
            species = species
    }

    call copy as copy_annotations{
        input:
            destination = destination + "/annotations",
            files =[
                   #haplo.out,
                   haplo.summary
                   ]
    }
}

task haplo {
    input {
        File vcf
        File? vcf_tbi
        String name = "haplo_output.json"
        String species = "homo_sapiens"
        Boolean database = false
        File fasta
        #Boolean offline = true
        File ensembl_cache

    }

    #TODO add haplo plugins http://www.ensembl.org/info/docs/tools/haplo/script/haplo_plugins.html
    #NOTE: the container has weird user permissions, using ln-s to avoid issues

    command {
        set -e
        ln -s ~{vcf} .
        ~{"ln -s "+ vcf_tbi  +" ."}
        ln -s ~{fasta} .
        mkdir .vep
        ln -s .vep /opt/vep/.vep
        haplo --verbose --input_file ~{basename(vcf)} -o ~{name} --species ~{species} --fasta ~{basename(fasta)} \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache}
    }
    runtime {
        docker: "quay.io/antonkulaga/vep"
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