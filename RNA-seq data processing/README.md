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

# 2. Check if a SV-associated gene is differentially expressed

We split accessions into SV-present and SV-absent groups and then applied one-sided Wilcoxon test on FPKM vectors of this gene to obtain P values. In each data source, only SV-gene pairs with the size of both SV-present and SV-absent groups no less than 5 were considered. P values for considered SV-gene pairs were corrected with Benjamini-Hochberg procedure (Benjamini and Hochberg, 1995). This was done through a custom R script. 

To identify SV-present and SV-absent groups for each SV-gene pair, we first prepared genotypes for deletions and insertions in 106 accessions obtained from the SV calling process: INS.genotype_gene_sv_pair.txt and DEL.genotype_gene_sv_pair.txt. Each line in the file shows information of a SV, and each column shows genotypes in '1/1' and '0/0'. 

The FPKM table for 73 accessions is DTA73MX-mBsnpREF-200425-expressed.csv. Each row is a gene and each column is an accession. 
```
Rscript SV-associated_gene_differentially_expressed_test.R
```

# 3. DEGs analyses

To identify differentially expressed genes (DEGs) between the collected RNA-seq data on three tissues (leaf, ear, and stem) for 6 maize and teosinte accessions, we applied the five methods to first find DEGs in each tissue with the false discovery rate(FDR) at 5%: cuffdiff (Ghosh and Chan, 2016), edgeR (v 3.28.1) (Robinson et al., 2010), Limma (v 3.42.2) (Smyth, 2004), SAMseq (v 3.0) (Li and Tibshirani, 2013) and DESeq2 (v 1.26.0) (Love et al., 2014). 

To run cufdiff analyses, the bam files were passed to cufflinks and the following commands were used to identify DEGs:
  ```
ADD lines here!
  ```
For edgeR, Limma, SAMseq, and DESeq2, we installed their R packages. In order to perform DEGs, we first quantify gene expression levels for each accession using HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/count.html) through the following command:
  
  ```
samtools view leaf_maize1.bam | htseq-count --stranded=no - Zm-Mo17-REFERENCE-CAU-1.0_genomic.gtf > leaf_maize1_htseq.txt 
```
Note that the gtf format of Mo17 genes Zm-Mo17-REFERENCE-CAU-1.0_genomic.gtf is like the following:
  ```
#gtf-version 2.2
#!genome-build Zm-Mo17-REFERENCE-CAU-1.0
#!genome-build-accession NCBI_Assembly:GCA_003185045.1
chr1	Genbank	gene	2759	7959	.	-	.	gene_id "Zm00014a_008050"; gbkey "Gene"; gene "MTERF2_2"; gene_biotype "protein_coding"; locus_tag "Zm00014a_008050"; partial "true"; 
chr1	Genbank	exon	7220	7959	.	-	.	gene_id "Zm00014a_008050"; transcript_id "gnl|WGS:NCVQ|Zm00014a_008050_T1"; gbkey "mRNA"; gene "MTERF2_2"; locus_tag "Zm00014a_008050"; orig_protein_id "gnl|WGS:NCVQ|Zm00014a_008050_P1"; orig_transcript_id "gnl|WGS:NCVQ|Zm00014a_008050_T1"; partial "true"; product "Transcription termination factor MTERF2, chloroplastic"; exon_number "1"; 
```

For each tissue, concatenate HTSeq results into a table, with each row represent a gene in the gtf file and each column represent the gene expression level of each accession, resulting in combine_htseq.txt. Run the following commands to obtain the DEGs in each tissue:
  ```
Rscript edger.R
Rscript deseq2.R
Rscript limma.R
Rscript samseq.R
```
DEGs supported by at least two of the five methods were considered as DEGs in each tissue.