# RNAseq-VC-nf
Nextflow variant calling pipeline for RNA-seq data using FreeBayes and/or GATK

## Requirements

### Call variants
- [Nextflow](https://www.nextflow.io/)
- [fastp](https://github.com/OpenGene/fastp)
- [STAR](https://github.com/alexdobin/STAR)
- [Picard](https://broadinstitute.github.io/picard/) (also packaged with GATK)
- [SAMtools](http://www.htslib.org/)
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [FreeBayes](https://github.com/ekg/freebayes)

### Filter variants
- [BCFtools](http://www.htslib.org/doc/bcftools.html)

### Annotate variants
- [BCFtools](http://www.htslib.org/doc/bcftools.html)
- [SnpEff](http://snpeff.sourceforge.net/)

## Usage
- Step 1: trim reads and two-pass alignment with STAR
    `nextflow run call_variants.nf [nextflow options] --fq [fastq dir]`

- Step 2: modify/merge BAMs, call variants
    `nextflow run call_variants.nf [nextflow options] --fq [fastq dir] --bam [bam dir]`

- Step 3: filter variants
    `nextflow run filter_variants.nf [nextflow options] --vcf [vcf dir]`

- Step 4: annotate variants
    `nextflow run annotate_variants.nf [nextflow options] --vcf [vcf dir] --filter [integer of the last filtering step]`
