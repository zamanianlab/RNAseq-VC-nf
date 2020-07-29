# Notes

## Daily log

### 2020/07/09

- Success of GATK best-practices for [RNA-seq](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) through all steps up to actual variant calling:
    - Two-pass alignment with STAR, utilizing SJs identified in the first pass to inform mapping in the second pass
    - Merging of uBAMs with aligned BAMs
    - Marking duplicates
    - Split reads that cross SJs
    - Recalibrate base quality scores with known *Ae. aegypti* LVP variants from Ensemble v46 (VCF from v47 is corrupted)

- Now planning/reading on the best approach for actual variant calling with Haplotype Caller
    - Can I use HaplotypeCaller in GVCF mode (as described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode)?)
        - Some evidence:
            - [Cattle](https://jasbsci.biomedcentral.com/articles/10.1186/s40104-019-0359-0)
            - [Barley](https://www.frontiersin.org/articles/10.3389/fpls.2019.00577/full#h3)
            - [Maize](https://www.nature.com/articles/s41598-019-56645-y#Sec2)
    - Restrict HaplotypeCaller to expressed intervals
        - `--intervals`
        - created `generate_intervals` process
    - Examples of hard-filtering:
        - [Trout](https://link.springer.com/article/10.1186/s12864-017-3992-z#Sec10)
            - Qual By Depth (QD) 2.0, FisherStrand (FS) 60.0: RMS Mapping Quality (MQ) 40.0, MAF > 0.05
        - [Barley](https://www.frontiersin.org/articles/10.3389/fpls.2019.00577/full#h3)
            - > 1 read depth and no neighbor polymorphisms around 60 bp
        - [Maize](https://www.nature.com/articles/s41598-019-56645-y#Sec2)
            - SNPs with a read depth of less than 5 were eliminated

- Can follow [this approach](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl) taken from the GATK-workflows repo

- If using the GCVF approach
    - add `-ERC GVCF` to the HaplotypeCaller arguments
    - use CombineGVCFs to combine by strain
    - use GenotypeGVCFs to genotype the strains

- Need to use the intervals generated in `generate_intervals` for BaseRecalibrator

- I think it's probably best to merge BAM files from the two strains immediately after MarkDuplicates
    - In which case the GCVF approach isn't necessary
    - decided to merge right before HaplotypeCaller

- Done through variant calling, now need to hard filter, do stats, and predict effects!

### 2020/07/16

- added arguments to be able to run from different stages without relying on `nextflow -resume`

### 2020/07/23

- HaploType caller is far too slow (3 days to finish ~10% of a single chromosome)
- add a process for freebayes
  - if we go this way, should we use the same pre-processing steps?
  - 'Following the guide in https://github.com/ekg/freebayes, we run Freebayes with default parameters after sorting and marking duplicates by GATK, then applied the "QUAL > 20" filter as recommended in tutorial. FreeBayes provides the useful options of ‘-C’ (minimal observations in a single sample, default value = 2) and ‘-F’ (minimal proportion of alternative reads, default value = 0.2) to filter the variants. When we used the parameters of “-C 0 –F 0”, the number of called SNVs became unreasonably large (52110 SNVs for Chr1 in average).' (https://doi.org/10.1186/s13059-019-1863-4)
- talk to CHTC
- transfer merged BAMs to a faster server at Mizzou

### 2020/07/26

- add FreeBayes process

### 2020/07/27

- VC by FreeBayes finished, now move to filtration

###2020/07/28

- move to new GitHub repo

1. Get unique variants
2. Run snpEff
3. Analyze

###2020/07/29

- some test filter commands:

  bcftools view -H --threads 8 Ref.vcf.gz | wc -l --- 1725083
  bcftools filter --threads 8 -e 'QUAL < 20' -Ou Ref.vcf.gz | bcftools view -H | wc -l --- 1172313


  bcftools view -H --threads 8 Sus.vcf.gz | wc -l --- 1797338
  bcftools filter --threads 8 -e 'QUAL < 20' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 1161406
  bcftools filter --threads 8 -e 'QUAL < 30' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 1012317
  bcftools filter --threads 8 -e 'QUAL < 40' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 801822
  bcftools filter --threads 8 -e 'QUAL < 50' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 720248
  bcftools filter --threads 8 -e 'INFO/DP < 50' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 423598
  bcftools filter --threads 8 -e 'INFO/DP < 10' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 823035
  bcftools filter --threads 8 -e 'MAF < 0.2' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 618693
  vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s Sus.vcf.gz | wc -l --- 60495

  # average depth at variant sites
  bcftools view -H Sus.vcf.gz | cut -f8 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' | gawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' --- 431.801
  # Heng Li filter 423 + 4*sqrt(432) < DP < 423 + 4*sqrt(432) https://arxiv.org/pdf/1404.0929.pdf
  bcftools filter --threads 8 -e 'INFO/DP > 506' -Ou Sus.vcf.gz | bcftools view -H | wc -l --- 1615365


- final command:
bcftools filter --threads 8 -e 'QUAL < 30' -Ou Sus.vcf.gz | \
  bcftools filter --threads 8 -e 'INFO/DP < 10' -Ou | \
  bcftools filter --threads 8 -e 'MAF < 0.2' -Ou | \
  bcftools filter --threads 8 -e 'INFO/DP > 506' -Oz -o Sus_1_filter.vcf.gz

bcftools index Sus_1_filter.vcf.gz
bcftools view -r 1:1-20000000 -Oz -o Sus_2_filter.vcf.gz Sus_1_filter.vcf.gz

bcftools filter --threads 8 -e 'QUAL < 30' -Ou Ref.vcf.gz | \
  bcftools filter --threads 8 -e 'INFO/DP < 10' -Ou | \
  bcftools filter --threads 8 -e 'MAF < 0.2' -Ou | \
  bcftools filter --threads 8 -e 'INFO/DP > 479' -Oz -o Ref_1_filter.vcf.gz

bcftools index Ref_1_filter.vcf.gz
bcftools view -r 1:1-20000000 -Oz -o Ref_2_filter.vcf.gz Ref_1_filter.vcf.gz
