
  cd PATH_TO_YOUR_SNPEFf/snpEff/
  java -jar ./snpEff.jar build -gff3 -v Mo17_CAU_1_0
  
  cd your_data_PATH
  perl get_modified_vcf.pl your_sample.heter0.05_miss_0.5.vcf  modified_SNP_VCF_file
  java -Xmx16g -jar PATH_TO_YOUR_SNPEFF/snpEff/snpEff.jar eff -c PATH_TO_YOUR_SNPEFF/snpEff/snpEff.config Mo17_CAU_1_0  modified_SNP_VCF_file > modified_SNP_VCF_annotation_file
    
  perl get_sys_and_no_syn_SNP.pl modified_SNP_VCF_annotation_file nonsynonymous_VCF synonymous_VCF
