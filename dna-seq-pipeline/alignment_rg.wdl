version development

workflow alignment {
   
    input {
        Array[File]+ reads
        File reference
        String name
        String destination
        Boolean rg_use_source
        File? rg_source 
        String ID
        String LB
        String PL        
        String PU        
        String SM         
        Boolean use_gencore
        Int compression# = 9
        Int align_threads# = 12
        Int sort_threads# = 12
        Int max_memory_gb# = 36
        Int coverage_sampling# = 1000
        String? gencore_quality 
    }

    call get_rg {
        input:
            rg_use_source = rg_use_source,
            rg_source = rg_source,   
            ID = ID,
            LB = LB,
            PL = PL,       
            PU = PU,
            SM = SM
    }


    call minimap2 {
        input:
            reads = reads,
            reference = reference,
            name = name,
            rg = get_rg.rg,
            threads = align_threads,
            max_memory = max_memory_gb
    }

    if (!use_gencore) {
        call sambamba_markdup {
            input:
                unsorted_bam = minimap2.bam,
                threads = sort_threads,
                compression = compression,
                max_memory = max_memory_gb
        }

    }

    if (use_gencore) {
        call sambamba_sort {
            input:
                unsorted_bam = minimap2.bam,
                threads = sort_threads,
                compression = 1,
                max_memory = max_memory_gb
        }

        call gencore{
            input:
                reference = reference,
                sorted_bam = sambamba_sort.bam,
                name = name,
                quality = gencore_quality,
                coverage_sampling = coverage_sampling,
                max_memory = max_memory_gb
        }

        call sambamba_sort as sambamba_sort_output {  #indexing fails due to gencore "WARNING: The output will be unordered!"
            input:
                unsorted_bam = gencore.bam,
                threads = sort_threads,
                compression = compression,
                max_memory = max_memory_gb
        }
    }
    
   
    call copy as copy_alignment {
        input:
            destination = destination,
            files = if (use_gencore) 
            then [select_first([sambamba_sort_output.bam, get_rg.e]),
                  select_first([sambamba_sort_output.bai, get_rg.e]),
                  select_first([gencore.html, get_rg.e]),
                  select_first([gencore.json, get_rg.e]),
                  select_first([sambamba_sort.flagstat, get_rg.e])] 
            else [select_first([sambamba_markdup.bam, get_rg.e]), 
                  select_first([sambamba_markdup.bai, get_rg.e]), 
                  select_first([sambamba_markdup.flagstat, get_rg.e])]
            
    }


    output {
       File bam =  copy_alignment.out[0]
       File bai = copy_alignment.out[1]
       Array[File] all = copy_alignment.out
    }
}

task get_rg {
    input {
        Boolean rg_use_source
        File? rg_source 
        String ID
        String LB
        String PL        
        String PU        
        String SM         
    }

    command {
        touch error.txt
        if [ "${rg_use_source}" != 'true' ]; then
          echo ~{"@RG\\\\\\\\tID:" + ID + "\\\\\\\\tLB:" + LB + "\\\\\\\\tPL:" + PL + "\\\\\\\\tPU:" + PU + "\\\\\\\\tSM:" + SM}
        else
          samtools view -H ~{rg_source} | grep '^@RG' | sed 's/'$'\t''/\\\\t/g'
        fi
    }

    runtime {
        docker_cpu: "1"
        docker: "quay.io/biocontainers/samtools@sha256:141120f19f849b79e05ae2fac981383988445c373b8b5db7f3dd221179af382b" #1.11--h6270b1f_0
    }

    output {
      File e = "error.txt" #terrible hack
      String rg = read_string(stdout())
    }
}


task minimap2 {
    input {
        Array[File] reads
        File reference
        String rg        
        String name
        Int threads
        Int max_memory
    }

    command {
        echo ~{rg}
        minimap2 -R ~{rg} -ax sr -t ~{threads} -2 ~{reference} ~{sep=' ' reads} | samtools view -bS - > ~{name}.bam
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker: "quay.io/comp-bio-aging/minimap2@sha256:69e9515a0cb5b5e9f47c3d0f95700d5064f0db04f86d49b6626b66a012daf0a5" #latest
        maxRetries: 2
      }

    output {
      File bam = name + ".bam"
    }
}

task sambamba_markdup {
    input {
        File unsorted_bam
        Int threads
        Int max_memory
        Int compression
        Int gb_per_thread = 3
    }

    String name = basename(unsorted_bam, ".bam")

    command {
       ln -s ~{unsorted_bam} ~{name + ".unsorted.bam"}
       sambamba markdup -t ~{threads} -p ~{name + ".unsorted.bam"} ~{basename(unsorted_bam)}
       sambamba sort -m ~{gb_per_thread}G -t ~{threads} -l ~{compression} -p ~{basename(unsorted_bam)}
       sambamba flagstat -t ~{threads} -p ~{name + ".sorted.bam"} > ~{name + ".sorted.bam.flagstat"}
    }

    runtime {
        docker: "quay.io/biocontainers/sambamba@sha256:ae92faef4c53a632b2120dfffa7b6dcfe5366a0647e61bbbd6188aedc89da4e8" #:0.8.0--h984e79f_0
        maxRetries: 1
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker_swap: "~{gb_per_thread * (threads+1) * 2}G"
      }

    output {
      File bam = name + ".sorted.bam"
      File bai = name + ".sorted.bam.bai"
      File flagstat = name + ".sorted.bam.flagstat"
    }
}


task sambamba_sort{
    input {
        File unsorted_bam
        Int threads
        Int max_memory
        Int compression
        Int gb_per_thread = 3
    }

    String name = basename(unsorted_bam, ".bam")

    command {
       ln -s ~{unsorted_bam} ~{basename(unsorted_bam)}
       sambamba sort -m ~{gb_per_thread}G -t ~{threads} -l ~{compression} -p ~{basename(unsorted_bam)}
       sambamba flagstat -t ~{threads} -p ~{name + ".sorted.bam"} > ~{name + ".sorted.bam.flagstat"}
    }

    runtime {
        docker: "quay.io/biocontainers/sambamba@sha256:ae92faef4c53a632b2120dfffa7b6dcfe5366a0647e61bbbd6188aedc89da4e8" #:0.8.0--h984e79f_0
        maxRetries: 1
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{threads+1}"
        docker_swap: "~{gb_per_thread * (threads+1) * 2}G"
      }

    output {
      File bam = name + ".sorted.bam"
      File bai = name + ".sorted.bam.bai"
      File flagstat = name + ".sorted.bam.flagstat"
    }
}




task gencore {
    input {
        File reference
        File sorted_bam
        String name
        Int max_memory
        Int supporting_reads = 1
        Float ratio_threshold = 0.8
        String? quality #"--high_qual"
        Int coverage_sampling# = 1000
    }
    command {
        gencore --coverage_sampling ~{coverage_sampling} --ratio_threshold=~{ratio_threshold} -s ~{supporting_reads} ~{quality} -i ~{sorted_bam} -o ~{name}.bam -r ~{reference}
        # samtools index ~{name}.bam  ~{name}.bam.bai
    }

    runtime {
        docker_memory: "~{max_memory}G"
        docker: "quay.io/comp-bio-aging/gencore@sha256:14b0da6c870766e04ea80a3d010ee593bf6a0bd071c5d4cdee002095a632a828"
    }

    output {
        File bam = name + ".bam"
        File html = "gencore.html"
        File json = "gencore.json"
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
