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

We mapped cis-eQTLs to scan for genetic markers that significantly affect gene expression in maize seedling. We selected common genetic markers (SVs, small InDels, and SNPs) with MAF >=0.05 for joint eQTL amapping. The joint eQTL analysis included 401,919 SVs, XX small InDels, and 7,153,302 SNPs. Gene expression levels of samples were quantile normalized, then subjected to inverse quantile normalization of each gene to control for outliers. We mapped cis-eQTLs with FastQTL v2.165 by using cis window of 1Mbp on either side of each gene. The covariates used in FastQTL included the first 3 genotyping principal components and 30 probabilistic estimation of expression residuals (PEER) factors. After running FastQTL, we performed gene-level multiple testing correction by using the Benjamini–Hochberg method at 5% FDR.


Create a genotype file: SV_SNP_INDEL.genotype.vcf.gz. Since only 73 accessions have RNA-seq data in seedling, we need to extract SVs, SNPs, and small InDels for these 73 samples and calculate genotyping PCAs. 

```
# The 73 samples were in wanted_sample

vcftools  --gzvcf GATK.108.snp.after_hard_filtered.vcf.gz  --keep wanted_sample   --recode --recode-INFO-all --stdout  > GATK.73.snp.after_hard_filtered.vcf.gz
vcftools  --gzvcf GATK.108.indel.after_hard_filtered.vcf.gz  --keep wanted_sample   --recode --recode-INFO-all --stdout  > GATK.73.indel.after_hard_filtered.vcf.gz

vcftools  --vcf DEL.SViper.out.1000.vcf.uniqe.nomissing.vcf  --keep wanted_sample   --recode --recode-INFO-all --stdout  > GATK.73.DEL.after_hard_filtered.vcf.gz
vcftools  --vcf INS.SViper.out.1000.vcf.uniqe.nomissing.vcf  --keep wanted_sample   --recode --recode-INFO-all --stdout  > GATK.73.INS.after_hard_filtered.vcf.gz

zcat GATK.73.DEL.after_hard_filtered.vcf.gz | perl hard_filter_GATK.pl - GATK.73.DEL.missrate05.het005.maf005.nobia.vcf.gz
zcat GATK.73.INS.after_hard_filtered.vcf.gz | perl hard_filter_GATK.pl - GATK.73.INS.missrate05.het005.maf005.nobia.vcf.gz

perl hard_filter_GATK.pl GATK.73.indel.after_hard_filtered.vcf.gz  GATK.73.indel.missrate05.het005.maf005.nobia.vcf.gz
perl hard_filter_GATK.pl GATK.73.snp.after_hard_filtered.vcf.gz    GATK.73.snp.missrate05.het005.maf005.nobia.vcf.gz

perl  convert_vcf_to_jvcf.pl  GATK.73.snp.missrate05.het005.maf005.nobia.vcf.gz   GATK.73.snp.missrate05.het005.maf005.nobia.jvcf.gz
perl  convert_vcf_to_jvcf_for_INDEL.pl GATK.73.indel.missrate05.het005.maf005.nobia.vcf.gz GATK.73.indel.missrate05.het005.maf005.nobia.jvcf.gz
perl  convert_vcf_to_jvcf_for_SV.pl GATK.73.INS.missrate05.het005.maf005.nobia.vcf.gz   GATK.73.INS.missrate05.het005.maf005.nobia.jvcf.gz
perl  convert_vcf_to_jvcf_for_SV.pl GATK.73.DEL.missrate05.het005.maf005.nobia.vcf.gz   GATK.73.DEL.missrate05.het005.maf005.nobia.jvcf.gz

perl convert_jvcf_into_ped_and_map_Huang2.pl GATK.73.snp.missrate05.het005.maf005.nobia.jvcf.gz GATK.73.snp.missrate05.het005.maf005.nobia.ped GATK.73.snp.missrate05.het005.maf005.nobia.map
perl convert_jvcf_into_ped_and_map_Huang2.pl GATK.73.indel.missrate05.het005.maf005.nobia.jvcf.gz GATK.73.indel.missrate05.het005.maf005.nobia.ped GATK.73.indel.missrate05.het005.maf005.nobia.map
perl convert_jvcf_into_ped_and_map_Huang2.pl GATK.73.INS.missrate05.het005.maf005.nobia.jvcf.gz   GATK.73.INS.missrate05.het005.maf005.nobia.ped GATK.73.INS.missrate05.het005.maf005.nobia.map
perl convert_jvcf_into_ped_and_map_Huang2.pl GATK.73.DEL.missrate05.het005.maf005.nobia.jvcf.gz   GATK.73.DEL.missrate05.het005.maf005.nobia.ped GATK.73.DEL.missrate05.het005.maf005.nobia.map

plink --vcf GATK.108.snp.hard_filter_missrate05.het005.maf005.nobia --allow-extra-chr --make-bed --out GATK.108.snp.hard_filter_missrate05.het005.maf005.nobia
plink --file GATK.73.DEL.missrate05.het005.maf005.nobia --allow-extra-chr --make-bed --out GATK.73.DEL.missrate05.het005.maf005.nobia 
plink --file GATK.73.INS.missrate05.het005.maf005.nobia --allow-extra-chr --make-bed --out GATK.73.INS.missrate05.het005.maf005.nobia
plink --file GATK.73.indel.missrate05.het005.maf005.nobia --allow-extra-chr --make-bed --out GATK.73.indel.missrate05.het005.maf005.nobia
plink --file GATK.73.snp.missrate05.het005.maf005.nobia   --allow-extra-chr --make-bed --out GATK.73.snp.missrate05.het005.maf005.nobia

DEL=GATK.73.DEL.missrate05.het005.maf005.nobia
INS=GATK.73.INS.missrate05.het005.maf005.nobia
SNP=GATK.73.snp.missrate05.het005.maf005.nobia
Indel=GATK.73.indel.missrate05.het005.maf005.nobia
awk '{print $2"\t"$2"_DEL"}' $DEL\.bim > ID_DEL.lst
awk '{print $2"\t"$2"_INS"}' $INS\.bim > ID_INS.lst
awk '{print $2"\t"$2"_Indel"}' $Indel\.bim > ID_Indel.lst

plink --bfile $DEL --update-map ID_DEL.lst --update-name --make-bed --out DEL --allow-extra-chr 
plink --bfile $INS --update-map ID_INS.lst --update-name --make-bed --out INS --allow-extra-chr  
plink --bfile $Indel --update-map ID_Indel.lst --update-name --make-bed --out Indel --allow-extra-chr 

INS=INS
DEL=DEL
SNP=GATK.73.snp.missrate05.het005.maf005.nobia
Indel=Indel

echo -e "INS.bed\tINS.bim\tINS.fam" > batch_list
echo -e "DEL.bed\tDEL.bim\tDEL.fam" >> batch_list
echo -e "Indel.bed\tIndel.bim\tIndel.fam" >> batch_list
echo -e "$SNP.bed\t$SNP.bim\t$SNP.fam" >> batch_list

plink --noweb  --merge-list batch_list --make-bed --out SNP_INDEL_SV  --allow-extra-chr
tabix -p vcf SV_SNP_INDEL.genotype.vcf.gz 
```

Prepare covariate.txt.gz file
```
# 3 genotyping PCAs

# Quantile normlization. The result file is seedling.normalized.expression.txt.

python quantile_normalization.py DTA73MX-mBsnpREF-200425-expressed.csv seedling.normalized --expression_threshold 0.1 --min_samples 10

# Calculate PEER factor. The results files are: seedling_PEER_covariates.txt, seedling_PEER_residuals.txt, and seedling_PEER_alpha.txt. 

Rscript run_peer.R

## Sort the columns in seedling_PEER_covariates.txt, make sure the sample order was the same in the genotype vcf.gz file.
sed '1d' seedling_PEER_covariates.txt | perl sort_phenotype_table.pl - seedling.normalized.expression.txt DEL.fam |cat header -  > seedling_PEER_covariates.txt_sorted

head -n4 PCA_covariants > covariates.txt
tail -n+2 seedling_PEER_covariates.txt_sorted | cat covariates.txt - | bgzip > covariates.txt.gz
```

Prepare a phenotype file: sorted.phenotype.bed.gz
```
# Sort the columns in seedling.normalized.expression.txt, make sure the sample order was the same in the genotype vcf.gz file. 

table=seedling.normalized.expression.txt
header=seedling.normalized.expression.txt
wanted_sort=DEL.fam

perl sort_phenotype_table.pl $table $header $wanted_sort > seedling.normalized.expression_sorted

# prepare correct phenotype files. This takes seedling.normalized.expression_sorted as input, and outputs phenotypes.bed. 

python add_gene_coordiates.py 

# index and zip phenotype files
# Remove genes that are not in the chr1---chr10. Input: phenotypes.bed,   output: phenotypes.bed.new. 
python remove_files.py  

sed -i s/chr//g phenotypes.bed.new 
head -1 phenotypes.bed.new > sorted.phenotype.bed
tail -n+2 phenotypes.bed.new | sort -k1,1 -k2,2n >> sorted.phenotype.bed
bgzip sorted.phenotype.bed && tabix -p bed sorted.phenotype.bed.gz
```
Run FastQTL the first time, perform gene-level multiple testing correction by using the Benjamini–Hochberg method at 5% FDR. 
```
marker=SV_SNP_INDEL
chunks=11
FastQTL/bin/fastQTL.static --vcf $marker.genotype.vcf.gz --bed sorted.phenotype.bed.gz --permute 1000 --out permutations_covariates --cov covariates.txt.gz --commands $chunks commands.covariates.$chunks.txt

# run commands in parallel
python generate_fastQTL_commands.py commands.covariates.$chunks.txt 11

# combine results from all chunks
cat permutations_covariates.* | gzip -c - > permutations_covariates.all.chunks.txt.gz

#  Selsect genes by multiple testing correction by using the Benjamini–Hochberg method at 5% FDR. The output file is: permutations_covariates.all.chunks.benjamini.txt, inclduing 12,513 eQTLs. 
Rscript BH_test_fdr0.05.R
```

# 3. Fine-mapping of causal variants at eQTLs

We used CAVIAR (https://github.com/fhormoz/caviar) to untangle linkage disequilibrium (LD) between genomic markers to predict a causal variation for eQTLs. For each eQTL, we selected the 100 most significant SNPs and small InDels, and 1 SV based on the FastQTL nominal P values. We estimated pairwise LD of the 101 markers for each gene and ran CAVIAR for each eQTL, by using the t statistics, and signed r values of LD among the 101 markers with a causal set size of 1. As a result, CAVIAR output the most significant variation for each eQTL.

```
#  Run fastQTL again for each eQTL with --map-full, display nominal P values for each variation for each gene.
## with --map 1050 option, and with covariates
FastQTL/bin/fastQTL.static --vcf SV_SNP_INDEL.genotype.vcf.gz --bed sorted.phenotype.bed.gz --map-full --out map_covariates --covcovariates.txt.gz --commands $chunks commands.covariates.$chunks.txt

python generate_fastQTL_commands.py commands.covariates.$chunks.txt 11

cat map_covariates.* > map_covariates.allchrs.combined.txt

# For each eGene, select ONE SV and 100 most significant SNPs and InDels
# CAVIAR requres zsocre, LD score for every gene
# zscore = -qnorm(p/2) if beta >0, zscore = -qnorm(1-p/2) if beta < 0
# LD score is calculated with plink, signed r value was used

# generate file combined.all.chunks.sh
python generate_CAVIAR_commands.py combined.all.chunks.benjamini.txt 100   

# split as 10 tasks and run in parallel
python generate_fastQTL_commands.py combined.all.chunks.sh 10  

cat combined.all.chunks.chunk*.txt.caviar.results > combined.all.chunks.CAVIAR.results.txt
python turn_1col_to_2cols_caviar_results.py  ## since some eGenes have two markers in the CAVIAR output.

grep -E "DEL|INS" $folder/$marker/combined.all.chunks.CAVIAR.results.2cols.txt > $folder/$marker/combined.all.chunks.CAVIAR.results.SV.txt  # 312 eGenes have a SV as causal variation.
```

