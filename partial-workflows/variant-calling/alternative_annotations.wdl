version 1.0
# -------------------------------------------------------------------------------------------------
# Package Name: SnpEff
# Task Summary: Annotate variants in a VCF
# Tool Name: eff/ann
# Documentation: http://snpeff.sourceforge.net/SnpEff_manual.html
# -------------------------------------------------------------------------------------------------


task SnpEff {
  input {
    File ? java
    File ? snpeff
    File config

    String filename_prefix
    File input_file
    File ? input_idx_file

    String dataDir
    String reference_version = "hg19"
    String ? userString

    Array[String] modules = []
    Float memory = 8
    Int cpu = 1

    String output_filename = filename_prefix + '.snpeff.vcf'
  }

  Int jvm_memory = round(memory)

  command {
    set -Eeuxo pipefail;

    for MODULE in ~{sep=' ' modules}; do
        module load $MODULE
    done;

    ~{default="java" java} \
      -Xmx~{jvm_memory}g \
      -jar ~{default="snpeff" snpeff} eff \
      ~{userString} \
      -c ~{config} \
      -dataDir ~{dataDir} \
      ~{reference_version} \
      ~{input_file} > ~{output_filename};
  }

  output {
    File vcf_file = "~{output_filename}"
  }

  runtime {
    memory: memory * 1.5 + " GB"
    cpu: cpu
  }

  parameter_meta {
    java: "Path to Java."
    snpeff: "SnpEff jar file."
    config: "Specify config file."
    input_file: "VCF file to annotate."
    input_idx_file: "VCF file index (.tbi)"
    dataDir: "Override data_dir parameter from config file."
    reference_version: "Version of genome to use."
    userString: "An optional parameter which allows the user to specify additions to the command line at run time."
    memory: "GB of RAM to use at runtime."
    cpu: "Number of CPUs to use at runtime."
  }

  meta {
    author: "Michael A. Gonzalez"
    email: "GonzalezMA@email.chop.edu"
    snpeff_version: "4.3q"
    version: "0.1.0"
  }
}
