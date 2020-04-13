version development

workflow index_reference {
input {
       File reference
       String? genome_assembly
       String? species
    }

    call create_reference_dict {
        input: reference=reference
    }

    call create_reference_fai {

    }

    output {

    }
}


task create_reference_dict {
    input {
        File reference
        String? assembly
        String? species
    }
    command {
        gatk CreateSequenceDictionary -R ~{reference} ~{"" + assembly}
    }
    runtime {
     docker: "quay.io/biocontainers/gatk4@sha256:4dcf52066fbee23130bf46b969adf85b55f6aec2c8871aba8d20b4b0a7115e1b"
    }

    output {

    }
}

