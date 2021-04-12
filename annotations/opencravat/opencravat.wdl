version development

workflow opencravat {
 input {
    File modules
    File jobs
    File logs
    File vcf
 }

    call annotate{
        input : modules = modules, jobs = jobs, logs = logs, vcf = vcf
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

        command {
            set -e
            ln -s ~{vcf} ~{basename(vcf)}
            export TMPDIR=/tmp
            oc run ~{basename(vcf)} -l hg38
            echo "annotation of the ~{vcf} finished!"
        }

        runtime {
            docker: "karchinlab/opencravat:latest"
            docker_volume1: modules_volume
            docker_volume2: jobs_volume
            docker_volume3: logs_volume
        }
        output {

        }
    }