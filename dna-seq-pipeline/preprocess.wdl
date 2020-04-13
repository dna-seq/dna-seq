version development

workflow preprocess {
    input {
       File reference
       String? genome_assembly
       String? species
    }


}


# Generate Base Quality Score Recalibration (BQSR) model
task recalibration_model {
    input {
               File bam
               File bai
               String report
               File reference
               File referenceDict
               File referenceFai
               String memory = "12G"
               String javaXmx = "4G"
    }

    String report_path = basename(bam, ".bam")+"_report"

    command {
        set -e
        mkdir -p "$(dirname ~{report_path})"
        gatk --java-options -Xmx~{javaXmx} \
        BaseRecalibrator \
        -R ~{reference} \
        -I ~{bam} \
        --use-original-qualities \
        -O ~{report_path} \
    }

    output {
        File out = report_path
    }


    runtime {
        docker: "quay.io/biocontainers/gatk4@sha256:4dcf52066fbee23130bf46b969adf85b55f6aec2c8871aba8d20b4b0a7115e1b" #4.1.6.0--py38_0
        memory: "16G"
    }
}

task apply_recalibration {
    input {
        File bam
        File bai
        File recalibration_report
        File reference
        File referenceDict
        File referenceFai
        String memory = "12G"
        String javaXmx = "4G"
    }

    String filename = basename(bam, ".bam")+"_recalibared"

    command {
        set -e
        mkdir -p "$(dirname ~{filename}_report)"
        gatk --java-options -Xmx~{javaXmx} \
        ApplyBQSR \
        --create-output-bam-md5 \
        --add-output-sam-program-record \
        -R ~{reference} \
        -I ~{bam} \
        -O ~{filename}.bam \
        --use-original-qualities \
        -bqsr ~{recalibration_report} \
        --static-quantized-quals 10 \
        --static-quantized-quals 20 \
        --static-quantized-quals 30 \
    }

    output {
        File recalibratedBam = filename+".bam"
        File recalibratedBamIndex =  filename+".bai"
        File recalibratedBamMd5 =  filename+".md5"
        File report = filename+"_report"
    }

    runtime {
        docker: "quay.io/biocontainers/gatk4@sha256:4dcf52066fbee23130bf46b969adf85b55f6aec2c8871aba8d20b4b0a7115e1b" #4.1.6.0--py38_0
        memory: "16G"
    }

}
