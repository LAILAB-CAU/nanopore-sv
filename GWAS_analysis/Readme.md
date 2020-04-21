基于表型和基因型数据进行GWAS分析

0.software
===
需要安装软件 GCTA,plink,emmax
```BASH
    conda install -c biobuilds gcta
    conda install -c biobuilds plink
    tar zxvf emmax-beta-07Mar2010.tar.gz
  ```
  
  1.genotype
  ===
  这里需要将vcf转化成特定格式的基因型文件，并得到表征血缘关系的的IBS.matrix与表征群体分层的PCA结果
  
  step1.转化vcf文件为map,和ped,bed文件，并过滤maf和missing rate
  ---
```BASH
    perl convert_vcf_to_hapmap.pl ./your_sample.vcf.gz your_sample.hapmap 
    gzip  your_sample.hapmap 
    perl convert_hapmap_to_jvcf.pl your_sample.hapmap.gz >  your_sample.jvcf
    perl get_MAF_0.05_missing_value_0.5_SNP.pl your_sample.vcf your_sample.maf.0.05,missing_value.0.5.jvcf
    perl convert_jvcf_into_ped_and_map_format_new.pl your_sample.maf.0.05,missing_value.0.5.jvcf your_sample.maf.0.05,missing_value.0.5.ped your_sample.maf.0.05,missing_value.0.5.map
    PATH_TO_PLINK/plink --file  your_sample.maf.0.05,missing_value.0.5 --map3 --missing-genotype 0 --make-bed --noweb --out  your_sample.maf.0.05,missing_value.0.5
 ```
step2.进行PCA分析
---
  
```BASH
    PATH_TO_GCTA/gcta64 --bfile your_sample.maf.0.05,missing_value.0.5 --autosome-num 10 --autosome --make-grm --out your_sample.maf.0.05,missing_value.0.5
    PATH_TO_GCTA/gcta64 --grm your_sample.maf.0.05,missing_value.0.5 --pca 20 --out  your_sample.maf.0.05,missing_value.0.5 
```
 
step3.得到tfam格式的基因型文件和IBS.matrix
---
```BASH
    PATH_TO_PLINK/plink  --file  your_sample.maf.0.05,missing_value.0.5  --map3 --noweb --missing-genotype 0 --recode12 --output-missing-genotype 0 --transpose --out  your_sample.maf.0.05,missing_value.0.5
    PTHA_TO_EMMAX/emmax-kin -v -h -s -d 10  your_sample.maf.0.05,missing_value.0.5
    perl kinship_add_inbred_name.pl your_sample.maf.0.05,missing_value.0.5.ped your_sample.maf.0.05,missing_value.0.5.hIBS.kinf > your_sample.maf.0.05,missing_value.0.5.kinship
 ```
step4:对PCA结果按照tfam文件中的自交系顺序排序
---
```BASH
     perl ./get_the_PCA_coverites_according_to_tfam_format.pl your_sample.maf.0.05,missing_value.0.5.eigenvec your_sample.maf.0.05,missing_value.0.5.tfam PCA.your_sample.maf.0.05,missing_value.0.5.tfam_order
```


2.phenotype
  ===
  对表型数据按照tfam文件的顺序排序，表型数据格式如下
  ```
    FAM1	B73	71.4081544
    FAM2	C002	75.93424471
    FAM3	C003	78.3163975
    FAM4	C004	75.93424471
    FAM5	C005	79.03104334
```
  排序脚本
```BASH
    perl sort_phenotype.pl phenotype_file PATH_TO_YOUR_genotype/your_sample.maf.0.05,missing_value.0.5.tfam phenotype_file_sorted
```

3.RUN EMMAX
  ===
```BASH
  PTHA_TO_EMMAX/emmax  -v -d 10 -t PATH_TO_YOUR_GENOTYPE/your_sample.maf.0.05,missing_value.0.5 -p PATH_TO_YOUR_PHENOTYPE/phenotype_file_sorted -c  PATH_TO_YOUR_GENOTYPE/PCA.your_sample.maf.0.05,missing_value.0.5.tfam_order  -k PATH_TO_YOUR_GENOTYPE/your_sample.maf.0.05,missing_value.0.5.hIBS.kinf  -o ./your_sample.GWAS_result
```

   

  
