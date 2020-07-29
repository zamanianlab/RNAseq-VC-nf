#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data=config.data_location // brc_location (btc) or data_location (other)
aux=config.aux_location
output=config.houtput_location
aedesgenome=config.aedesgenome_location

large_core=config.large_core
small_core=config.small_core

// Parameters

params.dir = null
if( !params.dir ) error "Missing dir parameter"
println "dir: $params.dir"

// flag for final stringtie_table_counts process (--stc)
params.stc = false

////////////////////////////////////////////////
// ** - Pull in fq files (paired); change file regular expression as needed
////////////////////////////////////////////////

Channel.fromFilePairs(data + "${params.dir}/*_R{1,2}.fq.gz", flat: true)
        .set { fq_pairs }

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.html'
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*.json'

   input:
       set val(id), file(forward), file(reverse) from fq_pairs

   output:
       set id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_fq_pairs
       set file("*.html"), file("*.json")  into trim_log

   """
       fastp -i $forward -I $reverse -o ${id}_R1.fq.gz -O ${id}_R2.fq.gz -y -l 50 -h ${id}.html -j ${id}.json
   """
}
trimmed_fq_pairs.set { trimmed_reads_hisat }



////////////////////////////////////////////////
// ** - LOAD in Aedes HiSat2 index and geneset files
////////////////////////////////////////////////
// aedesgenome points to /mnt/genomes/Other/Aedes_aegypti/

geneset_stringtie = file("${aedesgenome}/annotation/geneset_h.gtf.gz")
hs2_indices = Channel.fromPath("${aedesgenome}/Hisat2_indexes/*.ht2").buffer(size:8)




////////////////////////////////////////////////
// ** - HiSat2/Stringtie pipeline
////////////////////////////////////////////////

// alignment and stringtie combined
process hisat2_stringtie {

    publishDir "${output}/expression", mode: 'copy', pattern: '**/*'
    publishDir "${output}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        set val(id), file(forward), file(reverse) from trimmed_reads_hisat
        file("geneset.gtf.gz") from geneset_stringtie
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs
        file("${id}/*") into stringtie_exp
        file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
        hisat2 -p ${large_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        zcat geneset.gtf.gz > geneset.gtf
        stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
        rm *.gtf
    """
}

////////////////////////////////////////////////
// ** - STRINGTIE table counts & final commands that run on output dirs (run last with --stc flag)
// **  add later:    grep -Hn 'reads\|overall' ${output}/expression/*.hisat2_log.txt  | awk '{print $1}' | sed 's/.hisat2_log.txt//g' | sed 's/%//g' > Hisat2_stats.txt
//    publishDir "${output}/expression", mode: 'copy', pattern: 'Hisat2_stats.txt'
//        file("Hisat2_stats.txt") into hisat2_stats
////////////////////////////////////////////////

prepDE = file("${aux}/scripts/prepDE.py")
process stringtie_counts_final {

    echo true

    publishDir "${output}/counts", mode: 'copy', pattern: '*.csv'

    cpus small_core

    when:
      params.stc

    output:
        file ("gene_count_matrix.csv") into gene_count_matrix
        file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
        python ${prepDE} -i ${output}/expression -l 100 -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}
