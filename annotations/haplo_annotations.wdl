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
        String destination
    }

    call haplo{
        input: vcf = vcf,
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
        String name = "haplo_output.json"
        String species = "homo_sapiens"
        Boolean database = false
        File fasta
        #Boolean offline = true
        File ensembl_cache

    }

    #TODO add haplo plugins http://www.ensembl.org/info/docs/tools/haplo/script/haplo_plugins.html

    command {
        set -e
        haplo --verbose --input_file ~{vcf} -o ~{name} --species ~{species} --fasta ~{fasta} \
        ~{if(database) then "--database" else  "--cache"} --dir_cache ~{ensembl_cache}
    }
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