version development

# production configuration
import "https://raw.githubusercontent.com/dna-seq/dna-seq/master/dna-seq-pipeline/alignment.wdl" as alignment
import "https://raw.githubusercontent.com/dna-seq/dna-seq/master/dna-seq-pipeline/simple_variant_calling.wdl" as vc
import "https://raw.githubusercontent.com/dna-seq/dna-seq/master/dna-seq-pipeline/deep_variant.wdl" as deep

# debug local configuration (uncomment for debugging)
#import "alignment.wdl" as alignment
#import "simple_variant_calling.wdl" as vc
#import "deep_variant.wdl" as deep


workflow dna_seq_pipeline {
    input {
        Array[File]+ reads
        String destination
        Boolean is_paired = true
        File reference #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa
        File reference_fai #i.e. Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
        File? reference_index #used only if bwa-mem2 is an aligner
        Int align_threads# = 12
        Int sort_threads# = 12
        Int variant_calling_threads = 16
        Int coverage_sampling# = 1000
        Int max_memory_gb# = 42
        String name
        String? quality
        Boolean simple_variant_calling = false
        Boolean markdup = false
        Int compression = 9
        String sequence_aligner = "bwa-mem2"
    }


    call alignment.alignment as align{
        input:
          reads = reads,
          reference = reference,
          destination = destination,
          name = name,
          align_threads = align_threads,
          sort_threads = sort_threads,
          coverage_sampling = coverage_sampling,
          gencore_quality = quality,
          max_memory_gb = max_memory_gb,
          compression = compression,
          markdup = markdup,
          sequence_aligner = sequence_aligner,
          reference_index = reference_index
    }


    if(simple_variant_calling) {
        call vc.simple_variant_calling as simple_variant_calling{
            input:
                bam = align.out.bam,
                bai = align.out.bai,
                destination = destination + "/variants",
                referenceFasta = reference,
                referenceFai = reference_fai,
                threads = variant_calling_threads,
                max_memory = max_memory_gb,
                name = name
        }
        File simple_SNP = simple_variant_calling.results_SNP
        File simple_CV = simple_variant_calling.variants_SV
    }


    call deep.DeepVariant as deepvariant {
        input:
            bam = align.out.bam,
            bai = align.out.bai,
            name = name,
            reference = reference,
            reference_fai = reference_fai,
            threads = variant_calling_threads,
            destination = destination + "/variants/deepvariant",
            mode = "WGS" #[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    }



    output {
        File alignment = align.out.bam
        File deep_SNP = deepvariant.vcf
        File deep_g_SNP = deepvariant.gvcf
        File? results_SNP = simple_SNP
        File? results_SV =  simple_CV
    }


}
