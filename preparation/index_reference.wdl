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
        docker: "quay.io/biocontainers/gatk4@sha256:7b0b112b595861b140cbebdec5a0534bea9c40ef8bea4b3927fcea7ec53f5f57" #4.1.9.0--py39_0
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
        docker: "quay.io/biocontainers/samtools@sha256:141120f19f849b79e05ae2fac981383988445c373b8b5db7f3dd221179af382b" #1.11--h6270b1f_0
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
