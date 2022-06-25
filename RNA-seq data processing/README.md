In this study, we process RNA-seq data from three different sources. We generated RNA-seq data on seedlings for 73 out of the 106 accessions, collected unpublished RNA-seq data on seeds for 45 out of the 106 accessions, and collected RNA-seq data on three tissues (leaf, ear, and stem) for 6 maize and teosinte accessions (Lemmon et al., 2014). For each data source, we did FPKM calculation for each gene based on RNA-seq data, as well as identified differentially expressed genes (DEGs). 

# 1. FPKM calculation for each gene
Step 1. For the RNA-seq data on seedlings for 73 out of the 106 accessions, we used SNP-substituted Mo17 reference by including SNPs obtained through SNP calling pipeline. The following script takes C002 as an example, we extracted SNPs for each accession for substitution:
  ```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
gatk VariantFiltration \
-V GATK.106.snp.hard_filter.vcf.gz \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter" \
-O GATK_106.snps_filtered-hom.vcf.gz
gatk SelectVariants -V GATK_106.snps_filtered-hom.vcf.gz  \
--set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -select-type SNP -restrict-alleles-to BIALLELIC \
-O t106-hard3-snp.vcf.gz
gatk SelectVariants  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -restrict-alleles-to BIALLELIC  -R $refgenome -V t106-hard3-snp.vcf.gz -sn  C002   -O C002-hard2-snp.vcf.gz &
  ```

Step 2. Create SNP-substituted Mo17 reference:
  ```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
refgenome=$folder/Zm-Mo17-REFERENCE-CAU-1.0.fa
python mask_genome_variants.py --alternate-only  -i $refgenome -o C002--basedMo17.fa -v C002-hard2-snp.vcf
```

Step 3. RNA-seq reads trimming and filtering:
  ```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
fastp -w 24 --detect_adapter_for_pe    --cut_front --cut_tail     --in1 ../C2_L2_A008.R1.fastq.gz       --in2  ../C2_L2_A008.R2.fastq.gz    --out1   C2_L2_A008.R1.fastq.gz         --out2   C2_L2_A008.R2.fastq.gz           --report_title "C2_L2_A008"     --html  C2_L2_A008.html
```
Step 4. Reads mapping:
  
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
Step 5. FPKM values calculation:
  Download .gff3 file from https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/. 
```
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
cufflinks -p 6 -o ./cufflinks-C002-  -u -G Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 C002--basedMo17.us.bam
```

# 2. Common eQTL mapping

We mapped cis-eQTLs to scan for genetic markers that significantly affect gene expression in maize seedling. We selected common genetic markers (SVs, small InDels, and SNPs) with MAF >=0.05 for joint eQTL amapping. The joint eQTL analysis included 401,919 SVs, XX small InDels, and 7,153,302 SNPs. Gene expression levels of samples were quantile normalized, then subjected to inverse quantile normalization of each gene to control for outliers. We mapped cis-eQTLs with FastQTL v2.165 by using cis window of 1Mbp on either side of each gene. The covariates used in FastQTL included the first 3 genotyping principal components and 30 probabilistic estimation of expression residuals (PEER) factors. After running FastQTL, we performed gene-level multiple testing correction by using the Benjamini–Hochberg method at 10% FDR.

Prepare covariate.txt.gz file
```
# 3 genotyping PCAs

# Quantile normlization. The result file is seedling.normalized.expression.txt.

python quantile_normalization.py DTA73MX-mBsnpREF-200425-expressed.csv seedling.normalized --expression_threshold 0.1 --min_samples 10

# Calculate PEER factor. The results files are: seedling_PEER_covariates.txt, seedling_PEER_residuals.txt, and seedling_PEER_alpha.txt. 

Rscript run_peer.R

head -n4 ./SV_project/PCA_covariants > ./SV_SNP_INDEL/covariates.txt
tail -n+2 ./SV_project/seedling_PEER_covariates.txt_sorted | cat ./SV_SNP_INDEL/covariates.txt - | bgzip > ./SV_SNP_INDEL/covariates.txt.gz
```

Create a genotype file as FastQTL input: SV_SNP_INDEL.genotype.vcf.gz
```
tabix -p vcf SV_SNP_INDEL.genotype.vcf.gz
```

Prepare a phenotype file as FastQTL input: sorted.phenotype.bed.gz
```
# prepare correct phenotype files. This takes seedling.normalized.expression_sorted as input, and outputs phenotypes.bed. 
python add_gene_coordiates.py 

# # index and zip phenotype files
# Remove genes that are not in the chr1---chr10. Input: phenotypes.bed,   output: phenotypes.bed.new. 
python remove_files.py  

sed -i s/chr//g phenotypes.bed.new 
head -1 phenotypes.bed.new > sorted.phenotype.bed
tail -n+2 phenotypes.bed.new | sort -k1,1 -k2,2n >> sorted.phenotype.bed
bgzip sorted.phenotype.bed && tabix -p bed sorted.phenotype.bed.gz
```

Check if a SV-associated gene is differentially expressed

We split accessions into SV-present and SV-absent groups and then applied one-sided Wilcoxon test on FPKM vectors of this gene to obtain P values. In each data source, only SV-gene pairs with the size of both SV-present and SV-absent groups no less than 5 were considered. P values for considered SV-gene pairs were corrected with Benjamini-Hochberg procedure (Benjamini and Hochberg, 1995). This was done through a custom R script. 

To identify SV-present and SV-absent groups for each SV-gene pair, we first prepared genotypes for deletions and insertions in 106 accessions obtained from the SV calling process: INS.genotype_gene_sv_pair.txt and DEL.genotype_gene_sv_pair.txt. Each line in the file shows information of a SV, and each column shows genotypes in '1/1' and '0/0'. 

The FPKM table for 73 accessions is DTA73MX-mBsnpREF-200425-expressed.csv. Each row is a gene and each column is an accession. 
```
Rscript SV-associated_gene_differentially_expressed_test.R
```

