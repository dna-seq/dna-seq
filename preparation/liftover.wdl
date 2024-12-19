version development

#liftover pipeline
#using https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-quick-start.md as reference

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow DeepVariant {
    input {
        File input_vcf
        File chain_file
        File target_toplevel_fa
        String name
        String destination
        Int max_memory
    }

    call liftover {
        input:  input_vcf = input_vcf,
            chain_file = chain_file,
            target_toplevel_fa = target_toplevel_fa,
            name=name,
            max_memory = max_memory
    }
    call files.copy as copy_liftover {
        input: files = [liftover.vcf, liftover.altvcf, liftover.vcf_reject, liftover.altvcf_reject, liftover.liftover_log, liftover.liftover_alt_log], destination = destination
    }

    output {
        File vcf = copy_liftover.out[0]
        File altvcf = copy_liftover.out[1]
        File vcf_reject = copy_liftover.out[2]
        File altvcf_reject = copy_liftover.out[3]
        File liftover_log = copy_liftover.out[4]
        File liftover_alt_log = copy_liftover.out[5]
    }
}
task liftover{
    input {
        File input_vcf
        File chain_file
        File target_toplevel_fa
        String name
        Int max_memory
    }
    String output_vcf = name+".liftover.vcf"
    String output_altswap = name+".liftover.altswap.vcf"
    String output_reject = name+".liftover.reject.vcf"
    String output_alt_reject = name+".liftover.altswap.vcf"
    String output_log = name+".liftover.log"
    String output_alt_log = name+".liftover.altswap.log"

    command {
        ln -s ~{input_vcf} .
        ln -s ~{chain_file} .
        ln -s ~{target_toplevel_fa} .

        picard -Xmx~{max_memory}G -Xms128m  \
        LiftoverVcf \
        R=~{target_toplevel_fa} \
        I=~{input_vcf} \
        O=~{output_vcf} \
        CHAIN=~{chain_file} \
        REJECT=~{output_reject} \
        RECOVER_SWAPPED_REF_ALT=False > ~{output_log}

        picard -Xmx~{max_memory}G -Xms128m  \
        LiftoverVcf \
        R=~{target_toplevel_fa} \
        I=~{input_vcf} \
        O=~{output_altswap} \
        CHAIN=~{chain_file} \
        REJECT=~{output_alt_reject} \
        RECOVER_SWAPPED_REF_ALT=True > ~{output_alt_log}

        bgzip -l 9  ~{output_vcf}
        bgzip -l 9  ~{output_altswap}
        bgzip -l 9  ~{output_reject}
        bgzip -l 9  ~{output_alt_reject}
    }

    runtime {
        docker_memory: "~{max_memory+1}G"
        docker: "quay.io/biocontainers/picard@sha256:573ec3f38ab84c12a619eb31cf15cc1e5d709e50c083486cce04faf14ede91d4" #2.27.5--hdfd78af_0
        maxRetries: 2
    }
    output {
        File vcf = output_vcf+".gz"
        File altvcf = output_altswap+".gz"
        File vcf_reject = output_reject+".gz"
        File altvcf_reject = output_alt_reject+".gz"
        File liftover_log = output_log
        File liftover_alt_log = output_alt_log
    }
}