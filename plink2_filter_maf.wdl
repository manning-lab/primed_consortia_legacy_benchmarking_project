# plink2_filter_maf.
#
# This WDL workflow filters plink2 pgen files by MAF and combines the variant lists from multiple chromosomes into a single output file
#
# Required inputs:
#   - Array of pgen files
#   - Array of pvar files
#   - Array of psam files
#   - Float maf
#
# Optional inputs:
#   - Chromosome code (passed to "--output-chr")
#   - File with samples to keep
#   - File with variants to keep
#
# Output files:
#   - Filtered pgen/pvar/psam arrays for each chromosome
#   - combined_variant_IDs_MAF_gt_{maf}.txt
version 1.0

workflow plink2_filter_maf {
  input {
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    Float maf

    String output_chr_code = "chrM"
    File? samples_keep
    File? variant_IDs_keep
  }

  scatter (i in range(length(pgen_files))) {
    call plink2_filter_maf as plink2_filter_maf_chrom {
      input:
        pgen_file = pgen_files[i],
        pvar_file = pvar_files[i],
        psam_file = psam_files[i],
        maf = maf,
        output_chr_code = output_chr_code,
        samples_keep = samples_keep,
        variant_IDs_keep = variant_IDs_keep
    }
  }

  call combine_pvar_variant_ids {
    input:
      pvar_files = plink2_filter_maf_chrom.filtered_pvar,
      maf = maf
  }

  output {
    Array[File] filtered_pgen = plink2_filter_maf_chrom.filtered_pgen
    Array[File] filtered_pvar = plink2_filter_maf_chrom.filtered_pvar
    Array[File] filtered_psam = plink2_filter_maf_chrom.filtered_psam
    File combined_variant_list = combine_pvar_variant_ids.variant_ids_file
  }

  meta {
    author: "Alisa Manning"
    email: "amanning@broadinstitute.org"
  }
}

task plink2_filter_maf {
  input {
    File pgen_file
    File pvar_file
    File psam_file
    Float maf
    String output_chr_code
    File? samples_keep
    File? variant_IDs_keep
  }

  Int disk_size = ceil(3 * size(pgen_file, "GB")) + 10
  String bgen_prefix = basename(pgen_file, ".pgen")

  command <<<
    set -euo pipefail

    ln -s ~{pgen_file} ~{bgen_prefix}.pgen
    ln -s ~{psam_file} ~{bgen_prefix}.psam
    ln -s ~{pvar_file} ~{bgen_prefix}.pvar

    plink2 \
      --pfile ~{bgen_prefix} \
      ~{if defined(samples_keep) then "--keep " + samples_keep else ""} \
      ~{if defined(variant_IDs_keep) then "--extract " + variant_IDs_keep else ""} \
      --maf ~{maf} \
      --make-pgen \
      --output-chr ~{output_chr_code} \
      --out ~{bgen_prefix}_filtered_maf_~{maf}
  >>>

  output {
    File filtered_pgen = "~{bgen_prefix}_filtered_maf_~{maf}.pgen"
    File filtered_pvar = "~{bgen_prefix}_filtered_maf_~{maf}.pvar"
    File filtered_psam = "~{bgen_prefix}_filtered_maf_~{maf}.psam"
  }

  runtime {
    docker: "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
    disks: "local-disk ~{disk_size} SSD"
    memory: "16G"
  }
}

task combine_pvar_variant_ids {
  input {
    Array[File] pvar_files
    Float maf
  }

  Int disk_size = ceil(size(pvar_files, "GB")) + 10

  command <<<
    set -euo pipefail

    for pvar in ~{sep=' ' pvar_files}; do
      grep -v '^#' "$pvar" | awk '{print $3}'
    done > combined_variant_IDs_MAF_gt_~{maf}.txt
  >>>

  output {
    File variant_ids_file = "combined_variant_IDs_MAF_gt_~{maf}.txt"
  }

  runtime {
    docker: "ubuntu:24.04"
    disks: "local-disk ~{disk_size} SSD"
    memory: "4G"
  }
}
