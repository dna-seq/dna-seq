
# This task takes a BAM file and runs parliament2 on it.
task parlament2 {

  input
  {
    File bam
    File bai
    File reference_fai
    File reference_fasta
  }

  String name = basename(bam, ".bam")

  command {
      /home/dnanexus/parliament2.py --bam ~{bam} --bai ~{bai} --fai ~{reference_fai} -r ~{reference_fasta} \
      --breakdancer --breakseq --manta --cnvnator --lumpy \
      --delly_deletion --delly_insertion --delly_inversion --delly_duplication \
      --genotype --svviz
  }

  runtime {
    docker: "dnanexus/parliament2:latest"
  }

  output {
    File out = "output"
    #Array[File] vcfs = glob("output/*.vcf")
    #Array[File] sv_caller_results = glob("output/sv_caller_results/*")
    #Array[File] svtyped_vcfs = glob("output/svtyped_vcfs/*.vcf")
  }

}