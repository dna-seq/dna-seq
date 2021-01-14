version development

workflow bam_to_fastq {

    input {
        File bam #750018002018_WGZ.bam
        Int threads = 28
        String destination
    }

    call samtools_sort_by_name {
        input:   bam = bam, threads = threads
    }

    call bam2fastq{
        input: bam = samtools_sort_by_name.out
    }

    call copy {
        input: files = bam2fastq.out,
        destination = destination
    }

    output {
        Array[File] out = copy.out
    }
}

task samtools_sort_by_name {
    input {
        File bam
        Int threads
    }

    String name = basename(bam, ".bam")

    command {
       samtools sort -n ~{bam} --threads ~{threads} -o ~{name}_sorted.bam
    }

    runtime {
        docker: "quay.io/biocontainers/samtools@sha256:141120f19f849b79e05ae2fac981383988445c373b8b5db7f3dd221179af382b" #1.11--h6270b1f_0
        maxRetries: 2
    }

    output {
        File out = name + "_sorted.bam"
    }
}

task bam2fastq{
    input {
            File bam #750018002018_WGZ.bam
    }

    command {
        bedtools bamtofastq -i ~{bam} -fq ~{basename(bam, ".bam")}_1.fq -fq2 ~{basename(bam, ".bam")}_2.fq
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools@sha256:02e198f8f61329f9eafd1b9fc55828a31020b383403adec22079592b7d868006" #2.29.2--hc088bd4_0
    }

    output {
        Array[File] out = [basename(bam, ".bam")+"_1.fq", basename(bam, ".bam")+"_2.fq"]
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
