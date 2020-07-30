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
params.filter = null

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

// Manually made snpEff genome index for the VectorBase version of Aedes_aegypti:
//    mkdir /home/linuxbrew/.linuxbrew/Cellar/snpeff/4.3t/share/snpeff/data/Aedes_aegypti_VB
//    cp cp /mnt/genomes/Other/Aedes_aegypti/genome.fa /home/linuxbrew/.linuxbrew/Cellar/snpeff/4.3t/share/snpeff/data/Aedes_aegypti_VB/sequences.fa
//    cp cp /mnt/genomes/Other/Aedes_aegypti/annotation/geneset_h.gtf /home/linuxbrew/.linuxbrew/Cellar/snpeff/4.3t/share/snpeff/data/Aedes_aegypti_VB/genes.gtf
//    sed -i "s/>\(.*\) dna:.*/>\1/" sequences.fa
//    wget https://beta.vectorbase.org/common/downloads/release-47/AaegyptiLVP_AGWG/fasta/data/VectorBase-47_AaegyptiLVP_AGWG_AnnotatedCDSs.fasta
//    wget https://beta.vectorbase.org/common/downloads/release-47/AaegyptiLVP_AGWG/fasta/data/VectorBase-47_AaegyptiLVP_AGWG_AnnotatedProteins.fasta
//    cat VectorBase-47_AaegyptiLVP_AGWG_AnnotatedCDSs.fasta | perl -pe "s/>(.*?) .*/>\1/" > cds.fa
//    cat VectorBase-47_AaegyptiLVP_AGWG_AnnotatedProteins.fasta | perl -pe "s/>(.*?)-P(.) .*/>\1-R\2/" > protein.fa

// Removed a handful of small contigs that were causing build errors

// snpEff build -gtf22 -v Aedes_aegypti_VB

process variant_annotate {

      publishDir "${output}/annotations", mode: 'copy', pattern: '*_ann.vcf'
      publishDir "${output}/annotations", mode: 'copy', pattern: '*.html'
      publishDir "${output}/annotations", mode: 'copy', pattern: '*.txt'

      input:
          tuple val(id), file(vcf) from input_vcf

      when:
          params.vcf && params.filter

      """
         bcftools view --threads ${large_core} -Ov ${vcf} | snpEff -v Aedes_aegypti_VB - > ${id}_ann.vcf
      """
}
