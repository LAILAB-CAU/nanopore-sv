
# 1.Get sample SNP information
sampleINFOget.sh
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
gatk VariantFiltration \
    -V /public1/home/sc30797/ont_data/ngs_bams/gatk/population/GATK.106.snp.hard_filter.vcf.gz \
    --genotype-filter-expression "isHet == 1" \
    --genotype-filter-name "isHetFilter" \
    -O GATK_106.snps_filtered-hom.vcf.gz
gatk SelectVariants -V GATK_106.snps_filtered-hom.vcf.gz  \
    --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -select-type SNP -restrict-alleles-to BIALLELIC \
    -O t106-hard3-snp.vcf.gz
/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk SelectVariants  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -restrict-alleles-to BIALLELIC  -R ../Zm-Mo17-REFERENCE-CAU-1.0.fa -V t106-hard3-snp.vcf.gz -sn  C002   -O C002-hard2-snp.vcf.gz &
```

# 2.Genome SNP substitution 

snpmask.sh
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
python mask_genome_variants.py --alternate-only  -i Zm-Mo17-REFERENCE-CAU-1.0.fa -o C002--basedMo17.fa -v C002-hard2-snp.vcf
```



# 3.Reads trimming and filtering

fastp.sh
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
fastp -w 24 --detect_adapter_for_pe    --cut_front --cut_tail     --in1 ../C2_L2_A008.R1.fastq.gz       --in2  ../C2_L2_A008.R2.fastq.gz    --out1   C2_L2_A008.R1.fastq.gz         --out2   C2_L2_A008.R2.fastq.gz           --report_title "C2_L2_A008"     --html  C2_L2_A008.html
```


# 4.Reads mapping

hisat2.sh
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
hisat2-build   -p 24  C002--basedMo17.fa C002--basedMo17-hisat2
hisat2 -p 24 -q --dta-cufflinks -x ../C002--basedMo17-hisat2  -1 C2_L2_A008.R1.fastq.gz      -2 C2_L2_A008.R2.fastq.gz      -S C002--basedMo17.sam          
samtools view -@ 6 -bS  C002--basedMo17.sam -o C002--basedMo17.bam 
samtools view -q 20 -@ 6 -bS C002--basedMo17.bam -o C002--basedMo17.u.bam 
samtools sort -@ 6  C002--basedMo17.u.bam -o C002--basedMo17.us.bam
samtools index  -@ 6   C002--basedMo17.us.bam
```

# 5.FPKM values calculation

cufflinks.sh
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
cufflinks -p 6 -o ./cufflinks-C002-  -u -G Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 C002--basedMo17.us.bam
```

