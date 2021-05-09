version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow opencravat {
    input {
        File modules
        File jobs
        File logs
        File vcf
        String destination
        Array[String] report_formats = ["text", "tsv",  "csv", "rdata", "excel", "pandas"]
    }


    call annotate{
        input : modules = modules, jobs = jobs, logs = logs, vcf = vcf
    }


    call files.copy as copy_annotations
    {
        input:
            destination = destination,
            files = select_all([annotate.sqlite, annotate.err, annotate.log, annotate.status, annotate.extra, annotate.crv, annotate.crm, annotate.crm, annotate.crs])
    }

    call report {
        input: modules = modules, jobs = jobs, logs = logs, sqlite = copy_annotations.out[0], formats = report_formats
    }

    output {
        File sqlite = copy_annotations.out[0]
    }
}

task annotate {
    input {
        String modules
        String jobs
        String logs
        File vcf
    }

    String modules_volume = modules + ":/mnt/modules"
    String jobs_volume = jobs + ":/mnt/jobs"
    String logs_volume = logs + ":/mnt/logs"
    String name = basename(vcf)

    command {
        set -e
        ln -s ~{vcf} ~{name}
        export TMPDIR=/tmp
        oc run ~{name} -l hg38
        echo "annotation of the ~{vcf} finished!"
    }

    runtime {
        docker: "karchinlab/opencravat:latest"
        docker_volume1: modules_volume
        docker_volume2: jobs_volume
        docker_volume3: logs_volume
    }

    output {
        File err = name + ".err"
        File log = name + ".log"
        File sqlite = name + ".sqlite"

        File? status = name + ".status.json"
        File? extra = name + ".extra_vcf_info.var"
        File? crv = name + ".crv"
        File? crm = name + ".crm"
        File? crs = name + ".crs"
    }
}

task report {
    input {
        File sqlite
        String modules
        String jobs
        String logs
        Array[String] formats = ["text", "tsv",  "csv", "rdata", "excel", "pandas"]
    }

    String modules_volume = modules + ":/mnt/modules"
    String jobs_volume = jobs + ":/mnt/jobs"
    String logs_volume = logs + ":/mnt/logs"

    String name = basename(sqlite)

    command {
        set -e
        ln -s ~{sqlite} ~{name}
        export TMPDIR=/tmp
        oc report ~{name} -t ~{sep=" " formats}
    }

    runtime {
        docker: "karchinlab/opencravat:latest"
        docker_volume1: modules_volume
        docker_volume2: jobs_volume
        docker_volume3: logs_volume
    }
}