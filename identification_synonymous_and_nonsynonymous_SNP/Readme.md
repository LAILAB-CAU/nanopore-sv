基于SNP的位置和碱基变化判断SNP的效应

0.Software
===
  需要软件snpEff  
```BASH
  unzip snpEff_latest_core.zip
 ```

1:构建基于Mo17基因组和全基因组注释文件的基因组数据库
===

step1.编辑snpEff.config
---

```BASH
vim PATH_TO_YOUR_SNPEFf/snpEff/snpEff.config 
```

在'Third party databases'处修改

```BASH
# Third party databases
#-------------------------------------------------------------------------------

# Mo17 genome, version Mo17_CAU_1_0
      Mo17_CAU_1_0.genome : Mo17
 ```
 
 step2.重命名文件和文件夹
 ---
    在PATH_TO_YOUR_SNPEFf/snpEff/data下创建genomes和Mo17_CAU_1_0文件夹    
    将基因组序列文件命名为Mo17_CAU_1_0.fa，且放置在genomes文件夹内
    将全基因组注释文件命名为genes.gff，且放置在Mo17_CAU_1_0文件夹内
    
step3.构建snpEff数据库
---
  
```BASH
  cd PATH_TO_YOUR_SNPEFf/snpEff/
  java -jar ./snpEff.jar build -gff3 -v Mo17_CAU_1_0
```

2.注释SNP
===
基于硬过滤后的vcf文件(过滤杂合率和缺失率)


step1，保留vcf的前五列，并注释SNP的类型
---
```BASH
    perl get_modified_vcf.pl your_sample.heter0.05_miss_0.5.vcf  modified_SNP_VCF_file
    java -Xmx16g -jar PATH_TO_YOUR_SNPEFF/snpEff/snpEff.jar eff -c PATH_TO_YOUR_SNPEFF/snpEff/snpEff.config Mo17_CAU_1_0  modified_SNP_VCF_file > modified_SNP_VCF_annotation_file
```

step2.得到同义非同义突变的SNP位点,
---
```BASH
    perl get_sys_and_no_syn_SNP.pl modified_SNP_VCF_annotation_file nonsynonymous_VCF synonymous_VCF
 ```
 
 
