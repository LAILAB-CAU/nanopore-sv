# 0. Prepare Software
Install the latest_core version of SnpEff from https://sourceforge.net/projects/snpeff/files/. 
```
unzip snpEff_latest_core.zip
 ```
# 1. Build the reference genome and annotation file for Mo17
Step 1. Edit snpEff.config:
```
vim PATH_TO_YOUR_SNPEFf/snpEff/snpEff.config 
```
Change text at the 'Third party databases' as:
```
# Third party databases
#-------------------------------------------------------------------------------

# Mo17 genome, version Mo17_CAU_1_0
      Mo17_CAU_1_0.genome : Mo17
```

Step 2. Rename file and folder names:

Create /genomes and /Mo17_CAU_1_0 folders under PATH_TO_YOUR_SNPEFf/snpEff/data. Change the name of Mo17 reference genome (downloaded from https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/) into Mo17_CAU_1_0.fa and put it under /genomes folder. The Mo17 gene annotation file was renamed to genes.gff and was put under /Mo17_CAU_1_0 folder. 

Step 3. Build the snpEff database:
```
cd PATH_TO_YOUR_SNPEFf/snpEff/
java -jar ./snpEff.jar build -gff3 -v Mo17_CAU_1_0
```

# 2. Annotate SNPs derived from our SNP calling pipeline
Note that the SNPs used here was SNPs after the hard filtering step, only bi-allelic SNPs with missing rate < 50%, minor allele frequency (MAF) < 5%, and Heterozygosity rate < 5% were used.

Step 1. Keep the first five columns of your vcf file and annotate the type of SNPs:

```
perl get_modified_vcf.pl your_sample.heter0.05_miss_0.5.vcf  modified_SNP_VCF_file
java -Xmx16g -jar PATH_TO_YOUR_SNPEFF/snpEff/snpEff.jar eff -c PATH_TO_YOUR_SNPEFF/snpEff/snpEff.config Mo17_CAU_1_0  modified_SNP_VCF_file > modified_SNP_VCF_annotation_file
```
Step 2. Obtain the synonymous and non-synonymous SNPs:
```
perl get_sys_and_no_syn_SNP.pl modified_SNP_VCF_annotation_file nonsynonymous_VCF synonymous_VCF
 ```
 
 
