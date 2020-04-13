version development


workflow dna_seq_pipeline {
    input {
        Array[File]+ reads
        String destination
        #Boolean is_paired
        File genome
        String name
        Int threads = 3
    }

    call fastp as cleaning { input: reads = reads,  adapters = ["AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"], q = 32 }

    call copy as copy_quality {
        input:
            destination = destination + "/cleaned",
            files = [cleaning.reads_cleaned[0], cleaning.reads_cleaned[1], cleaning.report_json, cleaning.report_html]
    }

    call minimap2 {
        input:
            reads = cleaning.reads_cleaned,
            reference = genome,
            name = name,
            threads = 4
    }

    call samtools_conversion {
        input:
            sam = minimap2.out
    }

    call samtools_sort {
        input:
            bam = samtools_conversion.out,
            threads = 4
    }

    call gencore{
        input:
            reference = genome,
            sorted_bam = samtools_sort.out,
            name = name
    }

    call samtools_sort as sort_gencore{
        input:
            bam = gencore.out,
            threads = 4
    }

    call copy as copy_alignemnt {
        input:
              destination = destination + "/aligned",
              files = [sort_gencore.out, gencore.html, gencore.json]
    }


    output {
        File out =  sort_gencore.out
        File html = gencore.html
        File json = gencore.json
    }
}