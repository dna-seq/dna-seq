version development

workflow index_reference {
input {
       File reference
       String? assembly_name
       String? species
       String destination
    }

    call create_reference_dict {
        input: reference=reference, assembly = assembly_name, species = species
    }

    call create_reference_fai {
         input: reference = reference
    }

    call copy {
        input: destination = destination,
        files =[create_reference_dict.out, create_reference_fai.out]
    }

    output {
        File dict = copy.out[0]
        File fai = copy.out[1]
    }
}


task create_reference_dict {
    input {
        File reference
        String? assembly
        String? species
    }

    String name = sub(basename(reference, ".fasta"), ".fa", "")

    command {
        set -e
        gatk CreateSequenceDictionary -R ~{reference} -O ~{name}.dict ~{"--SPECIES " + species} ~{"--GENOME_ASSEMBLY " + assembly}
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4@sha256:4dcf52066fbee23130bf46b969adf85b55f6aec2c8871aba8d20b4b0a7115e1b"
    }

    output {
        File out = name + ".dict"
    }
}

task create_reference_fai {
    input {
        File reference
    }

    String name = basename(reference)

    command {
        ln -s ~{reference} ~{name}
        samtools faidx ~{name}
    }
    runtime {
        docker: "quay.io/biocontainers/samtools@sha256:97b9627711c16125fe1b57cf8745396064fd88ebeff6ab00cf6a68aeacecfcda" #1.2-0
    }

    output {
       File out = name+".fai"
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
