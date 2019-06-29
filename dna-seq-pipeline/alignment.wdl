version development

import "common.wdl" as common

workflow alignment {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
    }
}


task minimap2 {
    input {
        Array[File] reads
        File reference
        String name
        Int threads
    }

    command {
        minimap2 -ax sr  -t ~{threads} -2 ~{reference} ~{sep=' ' reads} > ~{name}.sam
    }

    runtime {
        docker: "docker pull quay.io/biocontainers/minimap2:2.17--h84994c4_0"
        maxRetries: 2
      }

    output {
      File out = name + ".sam"
    }
}

task samtools_conversion {
    input {
        File sam
    }

    String name = basename(sam, ".sam")

    command {
       samtools view -bS ~{sam} > ~{name}.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
        maxRetries: 2
      }

    output {
        File out = name + ".bam"
      }
}


task samtools_sort {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
       samtools sort ~{bam}  -o ~{name}_sorted.bam
    }

    runtime {
        docker: "biocontainers/samtools@sha256:6644f6b3bb8893c1b10939406bb9f9cda58da368100d8c767037558142631cf3"
        maxRetries: 2
      }

    output {
        File out = name + "_sorted.bam"
      }
}

task coverage {
    input {
        File bam
    }

    String name = basename(bam, ".bam")

    command {
        bedtools genomecov -bg -ibam ~{bam} > ~{name}.bedgraph
    }

     runtime {
            docker: "quay.io/biocontainers/bedtools@sha256:a0bb135afdec53be4b953a9a8efbc801cdb90706e6e63e11e3f60b06b8444f78" #2.23.0--he941832_1
            maxRetries: 2
          }

    output {
        File out = name + ".bedgraph"
    }
}


task picard_mark_duplicates {
    input {
        File bam
        String outputBamPath
        String metricsPath

        Int memory = 4
        Float memoryMultiplier = 3.0
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBamPath})
        picard -Xmx~{memory}G \
        MarkDuplicates \
        INPUT=~{bam} \
        OUTPUT=~{outputBamPath} \
        METRICS_FILE=~{metricsPath} \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CLEAR_DT="false" \
        CREATE_INDEX=true \
        ADD_PG_TAG_TO_READS=false \
        CREATE_MD5_FILE=true
    }

    output {
        IndexedBamFile out = object {
          file: outputBamPath,
          index: sub(outputBamPath, ".bam$", ".bai"),
          md5sum: outputBamPath + ".md5"
        }
        File metricsFile = metricsPath
    }

    runtime {
        docker: "quay.io/biocontainers/picard:sha256:b5750bf51f4223e0274430d07483c1be54e66fd5b96533c2d07a079c55da6972" #2.20.2--0
        memory: ceil(memory * memoryMultiplier)
    }
}