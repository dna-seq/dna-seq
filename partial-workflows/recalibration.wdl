version development

workflow recalibration {
    input {
        File reference
        File referenceDict
        File referenceFai
        File bam
        File bai
        String? genome_assembly
        String? species
    }

    call recalibration_model{
        input:
        bam = bam,
        bai = bai,
        reference = reference,
        referenceDict = referenceDict,
        referenceFai = referenceFai
    }


    call apply_recalibration{
        input:
        bam = bam,
        bai = bai,
        recalibration_report = recalibration_model.out,
        reference = reference,
        referenceDict = referenceDict,
        referenceFai = referenceFai
    }

    output {
        File recalibratedBam = apply_recalibration.recalibratedBam
        File recalibratedBamIndex =  apply_recalibration.recalibratedBamIndex
        File recalibratedBamMd5 =  apply_recalibration.recalibratedBamMd5
        File report = recalibration_model.out
    }
}


# Generate Base Quality Score Recalibration (BQSR) model
task recalibration_model {
    input {
               File bam
               File bai
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
        docker: "quay.io/biocontainers/gatk4@sha256:7b0b112b595861b140cbebdec5a0534bea9c40ef8bea4b3927fcea7ec53f5f57" #4.1.9.0--py39_0
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
    }

    runtime {
        docker: "a" #4.1.6.0--py38_0
        memory: "16G"
    }

}
