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
params.filter = necator_americanus_3_label

// Arguments allow you to choose which processes to run
// Note: when providing arguments, the vcf file parameter needs to include the
// full path in single-quotes
if( !params.vcf | !params.filter ) error "Missing parameter"


input_vcf = Channel.empty()

if (params.vcf) {
  input_vcf = Channel
                    .fromPath("${params.vcf}/*_${params.filter}_filter.bcf")
                    .map{file -> tuple(file.simpleName, file)}
} else exit 1, 'No path to data was provided.'

process variant_annotate {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_ann.vcf'

      input:
          tuple val(id), file(vcf) from input_vcf

      when:
          params.vcf && params.filter

      """
          snpEff download Aedes_aegypti

          snpEff -v -t ${large_core} -c /home/linuxbrew/.linuxbrew/opt/snpeff/share/snpeff/snpEff.config Aedes_aegypti ${vcf} > ${id}_ann.vcf
      """
}
