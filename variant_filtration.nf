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

input_vcf = Channel.empty()

if (params.vcf) {
  input_vcf = Channel
                    .fromPath("${params.vcf}/*.vcf.gz")
                    .map{file -> tuple(file.simpleName, file)}
} else exit 1, 'No path to data was provided.'

input_vcf.into { vcf1; vcf2}

// https://arxiv.org/pdf/1404.0929.pdf
// http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/
// http://www.ddocent.com/filtering/

process remove_indels {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(vcf) from vcf1

      output:
          tuple val(id), file("${id}_1_filter.bcf") into snp_filter

      when:
          params.vcf

      """
          bcftools view --threads 8 --exclude-types indels -Ob -o ${id}_1_filter.bcf ${vcf}
          bcftools index ${id}_1_filter.bcf
      """
}

process quality_filter {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(vcf) from snp_filter

      output:
          tuple val(id), file("${id}_2_filter.bcf") into quality_filter

      when:
          params.vcf

      """
          bcftools filter --threads 8 -e 'QUAL < 30' -Ob -o ${id}_2_filter.bcf ${vcf}
          bcftools index ${id}_2_filter.bcf
      """
}

process minimum_depth {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(vcf) from quality_filter

      output:
          tuple val(id), file("${id}_3_filter.bcf") into min_depth_filter

      when:
          params.vcf

      """
          bcftools filter --threads 8 -e 'INFO/DP < 10' -Ob -o ${id}_3_filter.bcf ${vcf}
          bcftools index ${id}_3_filter.bcf
      """
}

process qual_by_depth {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(og_vcf) from vcf2
          tuple val(id), file(filt_vcf) from min_depth_filter

      output:
          tuple val(id), file("${id}_4_filter.bcf") into qual_depth_filter

      when:
          params.vcf

      """
          AVG_DP=$(bcftools view -H ${og_vcf} | cut -f8 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' | gawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
          DP_THRESH=$(echo "\$AVG_DP + 4 * (sqrt(\$AVG_DP))" | bc)

          bcftools filter --threads 8 -e 'QUAL < \$DP_THRESH * 2 && INFO/DP > \$DP_THRESH' -Ob -o ${id}_4_filter.bcf ${filt_vcf}
          bcftools index ${id}_4_filter.bcf
      """
}

process allelic_balance {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_filter.bcf'

      input:
          tuple val(id), file(vcf) from qual_depth_filter

      output:
          tuple val(id), file("${id}_5_filter.bcf") into allele_filter

      when:
          params.vcf

      """
          bcftools filter --threads 8 -e 'MAF < 0.2' -Ob -o ${id}_5_filter.bcf ${vcf}
          bcftools index ${id}_5_filter.bcf
      """
}
