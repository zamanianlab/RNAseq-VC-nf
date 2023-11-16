#!/usr/bin/env nextflow

// Edit chtc.configuration!
input=params.input
output=params.output
aux=params.aux
genomes=params.genomes

large_core=config.large_core
small_core=config.small_core

big=params.big
small=params.small

// Parameters
params.aedesgenome = null
if( !params.aedesgenome ) error "Missing aedesgenome parameter"
println "aedesgenome: $params.aedesgenome"

params.dir = null
params.bam = null
params.vcf = null

// Arguments allow you to choose which processes to run
// Note: when providing arguments, the fq files only need to have the path after "/home/BIOTECH/zamanian/data/"
// and the bam files only need to have the path after "/home/BIOTECH/zamanian/GitHub/${repo}/h_output/"
// arguments need to be in single-quotes
if( !params.dir & !params.bam & !params.vcf ) error "Missing parameter"

input_dir = Channel.empty()
input_bam = Channel.empty()
input_vcf = Channel.empty()

if (params.dir && !params.bam) {
  input_dir = Channel
                    .fromFilePairs(input + "/${params.dir}/*_{1,2}.fq.gz", flat: true)
} else if (params.dir && params.bam) {
  input_bam = Channel
                    .fromPath(output + "/${params.bam}/*.bam")
                    .map{file -> tuple(file.baseName, file)}
  input_dir = Channel
                    .fromFilePairs(input + "/${params.dir}/*_{1,2}.fq.gz", flat: true)
} else if (params.vcf) {
  input_vcf = Channel
                    .fromPath(output + "/${params.vcf}/*.vcf.gz")
                    .map{file -> tuple(file.simpleName, file)}
} else exit 1, 'No path to data was provided.'

////////////////////////////////////////////////
// ** - Pull in fq files (paired); change file regular expression as needed
////////////////////////////////////////////////

 Channel.fromFilePairs(input + "/${params.dir}/*_{1,2}.fq.gz", flat: true)
        .set { fqs }

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trim_reads {

    publishDir "${output}/${params.dir}/trim_stats/", mode: 'copy', pattern: '*.{json,html}'

    cpus small
    tag { id }

    when:
        params.dir || params.bam

    input:
        tuple val(id), file(forward), file(reverse) from fqs

    output:
     //   tuple id, file("${id}_R1.fq.gz"), file("${id}_R2.fq.gz") into trimmed_reads_star1, trimmed_reads_star2, trimmed_reads_picard
        tuple id, file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into trimmed_reads_star1, trimmed_reads_star2, trimmed_reads_picard
        tuple file("*.html"), file("*.json") into trim_log

   """
       fastp -i $forward -I $reverse -w ${task.cpus} -o ${id}_1_trimmed.fq.gz -O ${id}_2_trimmed.fq.gz -y -l 150 -h ${id}.html -j ${id}.json
   """
}

////////////////////////////////////////////////
// ** - LOAD in Aedes STAR index and geneset files
////////////////////////////////////////////////
// aedesgenome points to /mnt/genomes/Other/Aedes_aegypti/

//star_index = Channel.fromPath(genomes + "${aedesgenome}/STARIndex/").collect() #genomes folder isn't set up quite yet

star_index = Channel.fromPath(input + "/${aedesgenome}/STAR_index/")
    .collect() 


////////////////////////////////////////////////
// ** - STAR PIPELINE
////////////////////////////////////////////////

// STAR parameters taken from https://doi.org/10.1038/s41586-018-0692-z
// First-pass alignment generates splice-junctions from all the samples
// (blind to strain/sample/replicate) to be used in the second-pass. We
// don't need BAMs here, just the SJs. NOTE: --genomeLoad can use LoadAndRemove
// here, but not in the second-pass.
process star_align_first {

    publishDir "${output}/${params.dir}/star", mode: 'copy', pattern: '*.Log.final.out'
    publishDir "${output}/${params.dir}/star", mode: 'copy', pattern: '*.flagstat.txt'
    publishDir "${output}/${params.dir}/counts", mode: 'copy', pattern: '*.ReadsPerGene.tab'
    publishDir "${output}/${params.dir}/sj", mode: 'copy', pattern: '*.tab'

    cpus big
    tag { id }
    maxForks 6

    when:
        params.dir && !params.bam

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_star1
        file index from star_index

    output:
        file("*out.tab") into mapped_splice_junctions
        tuple file("${id}.Log.final.out"), file("${id}.flagstat.txt") into alignment_logs_star
        file("${id}.ReadsPerGene.tab") into star_counts

    """
        STAR --readFilesIn ${forward} ${reverse} \
          --readFilesCommand zcat \
          --genomeDir ${index} \
          --genomeLoad LoadAndRemove \
          --runThreadN ${task.cpus} \
          --outFilterType BySJout \
          --alignIntronMax 1000000 \
          --alignMatesGapMax 1000000 \
          --outFilterMismatchNmax 999 \
          --outFilterMismatchNoverReadLmax 0.04 \
          --outFilterMultimapNmax 20 \
          --outFileNamePrefix ${id}. \
          --quantMode GeneCounts \
          --peOverlapNbasesMin 10 
          samtools flagstat ${id}.bam > ${id}.flagstat.txt
          cat ${id}.ReadsPerGene.out.tab | cut -f 1,2 > ${id}.ReadsPerGene.tab
    """
}

// The collectFile() function will concatenate all the junctions emited by sample
// in the mapped_splice_junctions channel. Will create a new queue channel called
// concatenated_junctions.

mapped_splice_junctions
   .collectFile(name: 'concatenated_junctions.txt', newLine: true)
   .first()
   .set { concatenated_junctions }

// Second-pass alignment uses the GTF (which was incorporated into the genome
// index) and the SJs generated from the first-pass. We emit aligned BAMs here.
// --outSAMattrRGline choices were based on this https://gatkforums.broadinstitute.org/gatk/discussion/6937/rna-seq-variant-calling-and-merging-of-sample-replicates
// SM: by strain (Ref or Sus)
// ID: by rep (Ref/Sus + Ctl/Inf/Naive + a/b/c)
// LB: by rep (Ref/Sus + Ctl/Inf/Naive + a/b/c)
// PL: ignore
process star_align_second {

    publishDir "${output}/star_log", mode: 'copy', pattern: '*Log*'
    publishDir "${output}/alignment_stats", mode: 'copy', pattern: '*alignment_stats.log'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_star2
        file index from star_index
        file junctions from concatenated_junctions

    output:
        file "${id}.Log.out" into alignment_logs
        tuple val(id), file("${id}.bam") into bam_files

    when:
        params.dir && !params.bam

    """
        SM=`echo ${id} | cut -c1-3 | tr -d '\n'`
        ID=${id}
        LB=${id}

        STAR \\
        --readFilesIn ${forward} ${reverse} \\
        --readFilesCommand zcat \\
        --genomeDir ${index} \\
        --sjdbFileChrStartEnd ${junctions} \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes NH HI AS NM MD \\
        --outSAMattrRGline ID:"\$ID" PL:illumina SM:"\$SM" LB:"\$LB" \\
        --outFilterType BySJout \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --outFilterMultimapNmax 20 \\
        --outFileNamePrefix ${id}.

        mv ${id}.Aligned.out.bam  ${id}.bam
        samtools flagstat -@ ${task.cpus} ${id}.bam > ${id}_alignment_stats.log
    """
}

////////////////////////////////////////////////
// ** - VARIANT CALLING PIPELINE
////////////////////////////////////////////////

// GATK likes to use both aligned and unaligned reads, so we first generate a
// uBAM from original FASTQs. NOTE: the rg variable will need to be adjusted
// based on the @RGs from the FASTQs (the regex could be made to be more robust.)
process picard_fastq_uBAM {

    input:
        tuple val(id), file(forward), file(reverse) from trimmed_reads_picard

    output:
        tuple id, file("${id}.ubam") into sorted_ubams

    when:
        params.bam

    """
        SM=`echo ${id} | cut -c1-3 | tr -d '\n'`
        ID=${id}
        LB=${id}

        gatk \
          --java-options \
          -Xmx6g \
          FastqToSam \
          -F1 ${forward} \
          -F2 ${reverse} \
          -O ${id}.ubam \
          -SM "\$SM" \
          -RG "\$ID" \
          -LB "\$LB" \
          -PL illumina \
          -SO queryname \
          -TMP_DIR ${output}
    """
}

// uBAMs and BAMs need to be QN sorted.
process picard_sort_bam {

    cpus large_core

    input:
        tuple val(id), file(bam) from input_bam

    output:
        tuple val(id), file("${id}_qn_sorted.bam") into sorted_bams

    when:
        params.bam

    """
        gatk \
          --java-options \
          -Xmx6g \
          SortSam \
          -I ${bam} \
          -O ${id}_qn_sorted.bam \
          -SO queryname \
          -TMP_DIR ${output}
    """
}

// This joins BAM and uBAM channels by ID.
joined_bams = sorted_bams.join(sorted_ubams)

// Merge BAM and uBAM.
process picard_merge {

    input:
        tuple val(id), file(bam), file(ubam) from joined_bams

    output:
        tuple val(id), file("${id}_merged.bam") into merged_bams

    when:
        params.bam

    """
        gatk \
          --java-options \
          -Xmx6g \
          MergeBamAlignment \
          -ALIGNED ${bam} \
          -UNMAPPED ${ubam} \
          -O "${id}_merged.bam" \
          -R "${aedesgenome}/genome.fa" \
          -TMP_DIR ${output}
    """
}

// Mark duplicates. NOTE: it's important that proper RG information is input
// into this process.
process picard_mark_duplicates {

    publishDir "${output}/picard_stats", mode: 'copy', pattern: '*_marked_dup_stats.txt'

    input:
        tuple val(id), file(bam) from merged_bams

    output:
        tuple val(id), file("${id}_duplicates.bam") into duplicate_bams
        file "${id}_marked_dup_stats.txt" into picard_logs

    when:
        params.bam

    """
        gatk \
          --java-options \
          -Xmx6g \
          MarkDuplicates \
          -I ${bam} \
          -O ${id}_duplicates.bam \
          -M ${id}_marked_dup_stats.txt \
          -TMP_DIR ${output}

    """
}

// GATK doesn't like reads that span SJs, so we split those into multiple reads.
process split_reads {

    input:
        tuple val(id), file(bam) from duplicate_bams

    output:
        tuple id, file("${id}_split.bam") into split_bams1, split_bams2

    when:
        params.bam

    """
        gatk \
          --java-options \
          -Xmx4g \
          SplitNCigarReads \
          -R "${aedesgenome}/genome.fa" \
          -I ${bam} \
          -O ${id}_split.bam \
          -TMP_DIR ${output}
    """
}

// Get a list of known Ae. aegypti LVP variants from Ensembl.
process fetch_variants {

    output:
        tuple file("known_variants.vcf.gz"), file("known_variants.vcf.gz.csi"), file("known_variants.vcf.gz.tbi") into known_variants

    when:
        params.bam

    """
        wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/variation/vcf/aedes_aegypti_lvpagwg/aedes_aegypti_lvpagwg.vcf.gz \
          -O known_variants.vcf.gz

        wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/variation/vcf/aedes_aegypti_lvpagwg/aedes_aegypti_lvpagwg.vcf.gz.csi \
          -O known_variants.vcf.gz.csi

        gatk IndexFeatureFile \
          -I known_variants.vcf.gz

    """
}

// Get transcript intervals from GTF (don't need to analyze untranscribed regions)
// and correct base-indexing when converting to BED. We have to use the shell
// block to allow BASH and Nextflow variables.
// Approach taken from https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
annotations = Channel.fromPath("${aedesgenome}/geneset.gff", type: 'file')

process generate_intervals {

      input:
          file(gff) from annotations

      output:
          file("transcript.interval_list") into (intervals1, intervals2, intervals3)
          file("genome_picard.dict") into (dict1, dict2, dict3)
          file("transcript_intervals.bed") into intervals_bed

      when:
          params.bam

      shell:
      //cat !{gff} | grep -E '\\stranscript\\s' | awk '{print $1 "\t" $4-1 "\t" $5}' > transcript_intervals.bed
      '''
          cat !{gff} | awk '$3 ~/^mRNA/' | awk '{print $1 "\t" $4-1 "\t" $5}' > transcript_intervals.bed

          picard CreateSequenceDictionary \
            R="!${aedesgenome}/genome.fa" \
            O=genome_picard.dict

          picard BedToIntervalList \
            I=transcript_intervals.bed \
            O=transcript.interval_list \
            SD=genome_picard.dict
      '''

}

// Recalibrate base quality scores. NOTE: again, proper RG data is important.
process recalibrate_bases {

    input:
        tuple val(id), file(bam) from split_bams1
        tuple file(vcf), file(index_csi), file(index_tbi) from known_variants
        file(intervals) from intervals1.first()
        file(reference_dict) from dict1.first()

    output:
        tuple id, file("${id}_recal_data.table") into brdt1, brdt2

    when:
        params.bam

    """
        gatk BaseRecalibrator \
          -R "${aedesgenome}/genome.fa" \
          --sequence-dictionary ${reference_dict} \
          -I ${bam} \
          -L ${intervals} \
          --known-sites ${vcf} \
          -O "${id}_recal_data.table"
    """
}

bam_recal_table = split_bams2.join(brdt1)

// Apply the recalibration.
process apply_recalibration {

    input:
        tuple val(id), file(bam), file(recal_table) from bam_recal_table
        file(intervals) from intervals2.first()
        file(reference_dict) from dict2.first()

    output:
        tuple stdout, file("${id}_recal.bam") into recal_bams

    when:
        params.bam

    """
        gatk ApplyBQSR \
          -R "${aedesgenome}/genome.fa" \
          --sequence-dictionary ${reference_dict} \
          -I ${bam} \
          -L ${intervals} \
          --bqsr-recal-file ${recal_table} \
          -O "${id}_recal.bam"

        echo ${id} | cut -c1-3 | tr -d '\n'
    """
}

// Generate summary plots to visualize improvements.
process plot_recalibration {

      publishDir "${output}/base_recalibration", mode: 'copy', pattern: '*_AnalyzeCovariates.pdf'

      input:
          tuple val(id), file(recal_table) from brdt2

      output:
          file("${id}_AnalyzeCovariates.pdf") into recal_plots

      when:
          params.bam

      """
          gatk AnalyzeCovariates \
            -bqsr ${recal_table} \
            -plots "${id}_AnalyzeCovariates.pdf"
 """
}

// In this case, we have many samples from only two strains, so it's best to
// merge all samples from like strain.
// -c combine RGs (don't add a discriminatory suffix)
// -p same as -c but for PGs
// samtools merge strips some RG info, so add generic RGs with Picard (SMs and RGs no longer matter)
process merge_bams_by_strain {

    input:
        tuple val(id), file(bams) from recal_bams.groupTuple()

    output:
        tuple id, file("${id}_rg_merged.bam") into merged_by_strain_bams_hc, merged_by_strain_bams_fb

    when:
        params.bam

    """
        samtools merge -cp ${id}_merged.bam ${bams.join(" ")}

        picard AddOrReplaceReadGroups \
          I=${id}_merged.bam \
          O=${id}_rg_merged.bam \
          RGID=${id} \
          RGLB=${id} \
          RGPL=ILLUMINA \
          RGPU=NA \
          RGSM=${id}

        samtools index ${id}_rg_merged.bam
    """
}

// Call variants.

process freebayes {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*.vcf'

      input:
          tuple val(id), file(bam) from merged_by_strain_bams_fb
          file(intervals) from intervals_bed.first()

      output:
          file("${id}.vcf") into vcfs_fb

      when:
          params.bam

      """
          ## for attempting freebayes-parallel
          # bedtools getfasta -fi "${aedesgenome}/genome.fa" -bed ${intervals} -fo transcript_intervals.fa
          # samtools faidx transcript_intervals.fa
          # cat transcript_intervals.bed | awk '{print \$1 ":" \$2 "-" \$3}' > transcript_intervals.txt

          freebayes \
            -f "${aedesgenome}/genome.fa" \
            -b ${bam} > "${id}.vcf"

          bgzip -@ 8 "${id}.vcf > ${id}.vcf.gz
      """

}


// process haplotype_caller {
//
//       publishDir "${output}/vcfs", mode: 'copy', pattern: '*.vcf.gz'
//       publishDir "${output}/vcfs", mode: 'copy', pattern: '*.vcf.gz.tbi'
//
//       input:
//           tuple val(id), file(bam) from merged_by_strain_bams_hc
//           file(intervals) from intervals3.first()
//           file(reference_dict) from dict3.first()
//
//       output:
//           file("${id}.vcf.gz") into vcfs
//           file("${id}.vcf.gz.tbi") into vcf_indices
//
//       when:
//           params.bam
//
//       """
//           gatk HaplotypeCaller \
//             -R "${aedesgenome}/genome.fa" \
//             -I ${bam} \
//             -L ${intervals} \
//             -O ${id}.vcf.gz \
//             --dont-use-soft-clipped-bases \
//             --standard-min-confidence-threshold-for-calling 20
//       """
// }

// Filter variants.

process filter_variants {

      publishDir "${output}/vcfs", mode: 'copy', pattern: '*_1_filter.vcf.gz'

      input:
          tuple val(id), file(vcf) from input_vcf

      output:
          file("${id}_1_filter.vcf.gz")

      when:
          params.vcf

      """
          bcftools filter --threads 8 -e 'QUAL < 30' -Ou ${vcf} | \
            bcftools filter --threads 8 -e 'INFO/DP < 10' -Ou | \
            bcftools filter --threads 8 -e 'MAF < 0.2' -Ou | \
            bcftools filter --threads 8 -e 'INFO/DP > 500' -Oz -o ${id}_1_filter.vcf.gz

          bcftools index ${id}_1_filter.vcf.gz

          bcftools view -r 1:1-20000000 -Oz -o ${id}_2_filter.vcf.gz ${id}_1_filter.vcf.gz
      """

}
