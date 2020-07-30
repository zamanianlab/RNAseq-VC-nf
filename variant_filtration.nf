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
// Note: when providing arguments, the fq files only need to have the path after "/home/BIOTECH/zamanian/data/"
// and the bam files only need to have the path after "/home/BIOTECH/zamanian/GitHub/${repo}/h_output/"
// arguments need to be in single-quotes
if( params.vcf ) error "Missing parameter"

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
          bcftools view -H --exclude-types indels -Ob -o ${id}_1_filter.bcf ${vcf}
          bcftools index ${id}_1_filter.bcf

          ### Quality filter
          # bcftools filter --threads 8 -e 'QUAL < 30' -Ob -o ${id}_2_filter.bcf ${id}_1_filter.bcf
          # bcftools index ${id}_2_filter.bcf

          ### Minimum depth filter
          # bcftools filter --threads 8 -e 'INFO/DP < 10' -Ob -o ${id}_3_filter.bcf ${id}_2_filter.bcf
          # bcftools index ${id}_3_filter.bcf

          ### Quality by depth filter
          # Get average depth
          # AVG_DP=$(bcftools view -H ${vcf} | cut -f8 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' | gawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
          # Calculate the average depth plus 4 times the square root of the average depth
          # DP_THRESH=$(echo "\$AVG_DP + 4 * (sqrt(\$AVG_DP))" | bc)
          # Filter
          # bcftools filter --threads 8 -e 'QUAL < \$DP_THRESH * 2 && INFO/DP > \$DP_THRESH' -Ob -o ${id}_4_filter.bcf ${id}_3_filter.bcf
          # bcftools index ${id}_4_filter.bcf

          # Allelic-balance filter
          # bcftools filter --threads 8 -e 'MAF < 0.2' -Ob -o ${id}_5_filter.bcf ${id}_4_filter.bcf
          # bcftools index ${id}_5_filter.bcf
      """

}
