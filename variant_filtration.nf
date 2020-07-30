#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data = config.data_location // brc_location (btc) or data_location (other)
aux = config.aux_location
output = config.output_location
aedesgenome = config.aedesgenome_location

large_core = config.large_core
small_core = config.small_core

// Parameters
params.vcf = null

// Arguments allow you to choose which processes to run
// Note: when providing arguments, the vcf file parameter needs to include the
// full path in single-quotes
if( !params.vcf ) error "Missing parameter"

input_fq = Channel.empty()
input_bam = Channel.empty()
input_vcf = Channel.empty()

if (params.vcf) {
  input_vcf = Channel
                    .fromPath("${params.vcf}/*.vcf.gz")
                    .map{file -> tuple(file.simpleName, file)}
} else exit 1, 'No path to data was provided.'


// https://arxiv.org/pdf/1404.0929.pdf
// http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/
// http://www.ddocent.com/filtering/

process filter_variants {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(vcf) from input_vcf

      output:
          tuple val(id), file("${id}_4_filter.bcf") into filtered_vcf

      when:
          params.vcf

      """
          ### Remove indels
          bcftools view --threads 8 --exclude-types indels -Ob -o ${id}_1_filter.bcf ${vcf}
          bcftools index ${id}_1_filter.bcf
      """

}
