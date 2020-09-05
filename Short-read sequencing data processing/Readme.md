# 0. Prepare softwares
To install softwares for short-read sequencing data processing, it is convenient to install anaconda python=3.7.4 first. Then install fastp, bwa, samtools, bcftools, GATK4, bgzip, tabix through conda:

```
conda install -c  bioconda fastp
conda install -c  bioconda bwa
conda install -c  bioconda samtools
conda install -c  bioconda bcftools
conda install -c  bioconda bgzip
conda install -c  bioconda tabix
```
Moreover, install [GATK4](https://github.com/broadinstitute/gatk/releases) beforehand.

# 1. Read preprocessing and mapping
Step 1. Prepare raw_fq.txt, the first column is accession number, and second/third column is the absolute path of the pair-end sequencing fastq file as follows:
```
B073 $folder/B73_1.fq.gz $folder/B73_2.fq.gz
B010 $folder/C010_1.fq.gz $folder/C010_2.fq.gz
...
```
Step 2. Use custom perl script to do quality control and generate batch scripts for BWA alignment

```
perl create_bash_file_for_bwa_index.pl raw_fq.txt
```

The generated batch script for BWA alignment looks as follows, take B073 as an example:
```
#!/bin/bash

# Specify the reference genome to use
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa  

bwa aln -t 32  $refgenome $folder/B073_fastp1.fq.gz -f B073_1.sai &
bwa aln -t 32  $refgenome  $folder/B073_fastp2.fq.gz -f B073_2.sai &
wait

bwa sampe -r "@RG\tID:B073\tPL:ILLUMINA\tLB:B073\tSM:B073"  $refgenome B073_1.sai B073_2.sai $folder/B073_fastp1.fq.gz $folder/B073_fastp2.fq.gz  | samtools view -b -@ 60  -q 20 -  | samtools sort -@ 60  - > B073.q20.sorted.bam
samtools index B073.q20.sorted.bam
```

# 2. Remove duplicates 
To mark duplicates in B073.q20.sorted.bam file, the gatk and samtools index commands were used. The output file is B073.q20.sorted.picard.bam.
```
gatk --java-options "-Xmx4G"  MarkDuplicates -I B073.q20.sorted.bam -O B073.q20.sorted.picard.bam  -M  B073.metrics  --VALIDATION_STRINGENCY SILENT 1>B073.mark 2>&1 && samtools index B073.q20.sorted.picard.bam &
wait
```

# 3. SNP calling

In step (1)-(4), we aimed to create a reference vcf file for BQSR process in step (5). 

## (1) Prepare the file list containing absolute path of *.q20.sorted.picard.bam files.
```
bamlist
```
## (2) Truncate the reference genome into 50 Mbp segments to parallize the SNP calling process. Write the start and end positions of the segments into Mo17_10chr_region.txt.sorted as follows:
```
chr10:1-50000000
chr10:50000001-100000000
chr10:100000001-149041351
```

## (3) Generate batch files for [bcftools](https://samtools.github.io/bcftools/bcftools.html) by using custom perl script:
```shell
perl create_calling_shell.pl Mo17_10chr_region.txt.sorted bcftools.sh
```
Each bcftools.sh files is like the following:
```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa  

bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:1-50000000  -f $refgenome -b bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_1-50000000.vcf.gz &
bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:50000001-100000000  -f $refgenome -b ./bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_50000001-100000000.vcf.gz &
bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:100000001-150000000  -f $refgenome -b ./bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_100000001-150000000.vcf.gz &
wait
```
## (4) concatenate SNPs in each genomic segment into whole-genome SNPs
In the last step, SNPs for each genomic segment ending up with *.vcf.gz, write down the absolute path of *.vcf.gz into Mo17_10chr_region_filelist.txt.sorted. Make sure the order of the files should follow the linear order of the reference genome. Generate index file for GATK4 as follows: 
```
bcftools concat -f Mo17_10chr_region_filelist.txt.sorted -o bcftools.vcf.gz -O z
tabix bcftools.vcf.gz
```
## (5) Base Quality Score Recalibration (BQSR) 
The introduction page for BQSR is https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-. Since there is no gold standard SNP set for maize, we use the vcf files called through bcftools in step (3) as the SNP training set. The BQSR procedure includes two parts: BaseRecalibrator and ApplyBQSR. 

Use custom perl script to generate batch scripts:
```
perl creat_bqsr_apply.pl
```
An example of the BaseRecalibrator script:
```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa

gatk BaseRecalibrator \
 -R $refgenome \
 -I B073.q20.sorted.picard.bam \
 --known-sites bcftools.vcf.gz \
 -O B073.sorted.markdup.recal_data.table && echo "** B073.sorted.markdup.recal_data.table done **" &
wait
```

An example of the ApplyBQSR script:
```shell
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
gatk ApplyBQSR \
 --bqsr-recal-file B073.sorted.markdup.recal_data.table \
 -R $refgenome \
 -I B073.q20.sorted.picard.bam \
 -O B073_marked.sorted.BQSR.bam && echo "** ApplyBQSR done **" &
wait
```

## (6) Calling gvcf with GATK4
GATK4 allows to call SNPs in a per-accession manner, and output results in *.g.vcf.gz files. Different from *.vcf.gz files, the genomic positions having the same genotypes with the reference genome will also be saved in the *.g.vcf.gz files. Similar to bcftools in step (3), we also parallelized gvcf calling process by truncating the reference genome into 106x10 = 1060 segments (the number of chromosome is 10 for maize). 


A custom perl script was used to generate commands to call gvcf:
```
perl creat_hc.pl
```

An example of the script to call gvcf
```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
gatk  --java-options "-Xmx4g" HaplotypeCaller \
  --emit-ref-confidence GVCF \
 -R $refgenome \
 -I B073_marked.sorted.BQSR.bam \
 -L chr1 \
 -O B073_marked.HC.chr1.g.vcf.gz && echo "** B073_marked.HC.1.g.vcf.gz done **" &
wait
```

## (7) Combine gvcf files and joint calling
Combine the gvcf results for 1060 genomic segments into *.g.vcf.gz for each chromosome, GenotypeGVCFs option was then used to obtain the vcf file as *.vcf.gz for each chromosome.

A custom perl script was used to generate batch combining script:
```
perl create_bash_file_for_Combine.pl
```

An example for the combining script:
```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
gatk  --java-options "-Xmx128g" CombineGVCFs  \
 -R $refgenome \
 -V B73_marked.HC.chr1.g.vcf.gz \
 -V C002_marked.HC.chr1.g.vcf.gz \
 -V C003_marked.HC.chr1.g.vcf.gz \
   ...
 -V Mo17_marked.HC.chr1.g.vcf.gz \
 -O $outdir/GATK.106.HC.chr1.g.vcf.gz && echo "** GATK.106.HC.chr1.g.vcf.gz done **" 
gatk  --java-options "-Xmx128g" GenotypeGVCFs  \
 -R $refgenome \
 -V $outdir/GATK.106.HC.chr1.g.vcf.gz \
 -O $outdir/GATK.106.HC.chr1.vcf.gz && echo "** GATK.106.HC.chr1.vcf.gz done **" 
```

## (8) Merge vcf files
Combine *.vcf.gz for each chromosome into a genome-wide vcf file. 

```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
$gatk MergeVcfs \
 -I GATK.106.HC.chr1.vcf.gz \
 -I GATK.106.HC.chr10.vcf.gz \
 -I GATK.106.HC.chr2.vcf.gz \
 -I GATK.106.HC.chr3.vcf.gz \
 -I GATK.106.HC.chr4.vcf.gz \
 -I GATK.106.HC.chr5.vcf.gz \
 -I GATK.106.HC.chr6.vcf.gz \
 -I GATK.106.HC.chr7.vcf.gz \
 -I GATK.106.HC.chr8.vcf.gz \
 -I GATK.106.HC.chr9.vcf.gz \
 -O GATK.106.HC.vcf.gz && echo "** MergeVcfs done **"
```

## (9) Filter the genome-wide SNPs
Step 1. Extract SNP info from the vcf file:
```
gatk SelectVariants -V GATK.106.HC.vcf.gz -select-type SNP -O GATK.106.Wraw-snp.vcf.gz
```

Step 2. [Filter variants either with VQSR or by hard-filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering):
```
gatk VariantFiltration \
    -V GATK.106.Wraw-snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O GATK.106.snps_filtered-2.vcf.gz
```
Step 3. Keep homogenous sites only:
```
gatk VariantFiltration \
    -V GATK.106.snps_filtered-2.vcf.gz \
    --genotype-filter-expression "isHet == 1" \
    --genotype-filter-name "isHetFilter" \
    -O GATK.106.snps_filtered-hom.vcf.gz

gatk SelectVariants -V GATK.106.snps_filtered-hom.vcf.gz  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -select-type SNP -restrict-alleles-to BIALLELIC -O t106-hard2-snp.vcf.gz
```

Step 4. Extract the SNPs for a certain accession (take B073 as an example):
```
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
gatk SelectVariants  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -restrict-alleles-to BIALLELIC  -R $refgenome -V t106-hard2-snp.vcf.gz -sn B073 -O B073-hard2.vcf 
```
