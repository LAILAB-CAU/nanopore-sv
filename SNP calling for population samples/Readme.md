[TOC]
# 0.software
需要fastp、bwa、samtools、bcftools、GATK4、bgzip、tabix

依赖conda安装，annaconda python=3.7.4
```shell
conda install -c  bioconda fastp
conda install -c  bioconda bwa
conda install -c  bioconda samtools
conda install -c  bioconda bcftools
conda install -c  bioconda bgzip
conda install -c  bioconda tabix
```
安装[GATK4](https://github.com/broadinstitute/gatk/releases)
```shell
unzip gatk-4.1.0.0.zip
```

# 1.质控与比对
首先，使用fastp对从测序公司返回的重测序fastq文件进行质控，并产生通过质控过滤的fastq文件，以及质控的jsonhehtml格式的报告。之后，使用bwa aln将经过质控后的fastq文件比对到参考基因组上，生成.sai格式的文件，再通过bwa sampe将成对的.sai转化为比对的结果文件，默认输出是sam文件，为了节省存储空间可以通过管道命令接samtools直接生成.bam文件。

## 1.整理出sample的编号及其双端的fastq文件路径
>raw_fq.txt示例：第一列是样品编号，第二三列分别是该样品双端测序文件的路径 
>>>B073	/public1/home/sc30797/ont_data/ngs/F19FTSNCKF1832_li201910141623401653RE0101/Clean/B73/B73_1.fq.gz	/public1/home/sc30797/ont_data/ngs/F19FTSNCKF1832_li201910141623401653RE0101/Clean/B73/B73_2.fq.gz<br/>
C010	/public1/home/sc30797/ont_data/ngs/F19FTSNCKF1832_li201910141623401653RE0101/Clean/C10/C10_1.fq.gz	/public1/home/sc30797/ont_data/ngs/F19FTSNCKF1832_li201910141623401653RE0101/Clean/C10/C10_2.fq.gz

## 2.使用perl脚本批量生成质控+bwa比对的脚本
```shell
perl create_bash_file_for_bwa_index.pl raw_fq.txt
```

生成的脚本示例
```shell
#!/bin/bash
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa

#'-t' 规定了运行本条命令使用的线程数，可以根据服务器情况调整；'-f'规定了生成的.sai文件的命名；命令最后的'&'表示将本条命令挂在后台,此时可以在前台运行其他命令
bwa aln -t 32  $refgenome  /public1/home/sc30797/wangzijian/sv/Supplement/C003_fastp1.fq.gz -f C003_1.sai &
bwa aln -t 32  $refgenome  /public1/home/sc30797/wangzijian/sv/Supplement/C003_fastp2.fq.gz -f C003_2.sai &
#'wait'表示要等待之前挂在后台的命令完成后再运行下一步的命令，我们必须要等双端的fastq文件比对结束才能进行sampe的运行
wait
# -r "@RG\tID:C003\tPL:ILLUMINA\tLB:C003\tSM:C003" 参数决定了比对文件的描述
#后续四个文件分别为sai文件1 sai文件2 经fastp质控后的fastq文件1 经fastp质控后的fastq文件2
#'|'管道代表将上一步的输出作为这一步的输入，'-'就代表的是这个输入
#'-b'指的是输出成.bam格式；'-q 20'指的是只保留比对质量大于20的比对结果，为保证比对质量一般的分析我们的比对结果都应该过滤q20；'-@ 60'指的是调用的线程数
bwa sampe -r "@RG\tID:C003\tPL:ILLUMINA\tLB:C003\tSM:C003"  $refgenome C003_1.sai C003_2.sai /public1/home/sc30797/wangzijian/sv/Supplement/C003_fastp1.fq.gz /public1/home/sc30797/wangzijian/sv/Supplement/C003_fastp2.fq.gz  | samtools view -b -@ 60  -q 20 -  | samtools sort -@ 60  - > C003.q20.sorted.bam
samtools index C003.q20.sorted.bam
```

# 2.MarkDuplicates去冗余
```shell
#上一步的.q20.sorted.bam比对文件去冗余，得到.q20.sorted.picard.bam文件。
#注意，因为bwa aln+sampe工具比较老，所以加入'--VALIDATION_STRINGENCY SILENT'参数防止程序不正常终止
#'&&'表明gatk和samtools index这两个命令是衔接的，gatk完成后即进行samtools index
#可以将多个材料的命令写入同一个shell脚本,使用'&'后台挂起的方式同步运行，注意计算节点的内存和线程数的上限
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
$gatk --java-options "-Xmx4G"  MarkDuplicates -I ./C003.q20.sorted.bam -O C003.q20.sorted.picard.bam  -M  C003.metrics  --VALIDATION_STRINGENCY SILENT 1>C003.mark 2>&1 && samtools index C003.q20.sorted.picard.bam &
wait
```

# 3.bcftools call snp
这个主要是为了GATK4的BQSR过程提供一个参考的vcf。此前这个工具在samtools中，新版本将其整合到bcftools里了。
## 1.准备.q20.sorted.picard.bam文件的路径列表
bamlist
## 2.准备参考基因组的分段文件
因为全基因组call snp耗时很长，将参考基因组分成片段
Mo17_10chr_region.txt.sorted
>Mo17示例（chr10）：
>>chr10:1-50000000<br/>
chr10:50000001-100000000<br/>
chr10:100000001-149041351


## 3.批量生成bcftools命令脚本
```shell
perl create_calling_shell.pl Mo17_10chr_region.txt.sorted bcftools.sh
```

[bcftools使用说明](https://samtools.github.io/bcftools/bcftools.html)

bcftools.sh示例
```shell
bcftools=~/wangzijian/software/bcftools-1.10.2/my_bcftools/bin/bcftools
#conda bcftools运行时总是报错缺lib文件，所以我先规定一下这个总是缺的lib文件，把它写到环境里
export LD_LIBRARY_PATH=/public1/home/sc30797/ont_data/ngs_bams/bcf_snp/lib:$LD_LIBRARY_PATH
#这里依然是'&'后台并行
$bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:1-50000000  -f /public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa -b ./bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_1-50000000.vcf.gz &
$bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:50000001-100000000  -f /public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa -b ./bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_50000001-100000000.vcf.gz &
$bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r  chr1:100000001-150000000  -f /public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa -b ./bamlist | bcftools call -cv -Oz --threads 21 -o  bcftools_106__chr1_100000001-150000000.vcf.gz &
wait
```

## 4.将call snp结束后的片段接成全基因组snp
准备所有片段call snp的结果文件.vcf.gz,将路径写在一个文件里***Mo17_10chr_region_filelist.txt.sorted***，一行一个，注意片段的顺序要严格按照基因组的顺序排列
```shell
bcftools=~/wangzijian/software/bcftools-1.10.2/my_bcftools/bin/bcftools
export LD_LIBRARY_PATH=/public1/home/sc30797/ont_data/ngs_bams/bcf_snp/lib:$LD_LIBRARY_PATH
$bcftools concat -f Mo17_10chr_region_filelist.txt.sorted -o bcftools.vcf.gz -O z
#这里会生成bcftools.vcf.gz的索引文件bcftools.vcf.gz.tbi,因为GATK4要读取vcf的话一般是要求有这样一个文件的
tabix bcftools.vcf.gz
```

# 4.BQSR
进行BQSR是为了校正碱基质量值，以去除测序过程中由各种情况造成的错误。BQSR过程将所有与参考基因组不一致的变异视为错配，但是会将提供的标准SNP位点的错配视为合法错配，没有变异的位点的错配视为真正的测序错误，进行训练。人类有已知的标准SNP数据库，玉米里没有这样的数据，一般处理这样的物种时，可以使用未经BQSR的bam文件call一版，用这一版vcf作为训练用的变异数据库，也可以在此基础上多次迭代，当然也可以使用不同call snp的工具取交集。我们用之前bcftools call的snp作为训练的数据库。

分为两步：BaseRecalibrator和ApplyBQSR

之前bcftools时准备好的.q20.sorted.picard.bam文件的路径列表也可以拿来用

perl脚本生成批量的BQSR命令
```shell
perl creat_bqsr_apply.pl
```

BaseRecalibrator脚本示例
```shell
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa
path=~/ont_data/ngs_bams/bqsr
$gatk BaseRecalibrator \
 -R $refgenome \
 -I /public1/home/sc30797/ont_data/ngs_bams/aln_bam2/C113.q20.sorted.picard.bam \
 --known-sites /public1/home/sc30797/ont_data/ngs_bams/bcf_snp/bcftools.vcf.gz \
 -O $path/C113.sorted.markdup.recal_data.table && echo "** C113.sorted.markdup.recal_data.table done **" &
wait
```

ApplyBQSR脚本示例
```shell
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa
path=~/ont_data/ngs_bams/bqsr
$gatk ApplyBQSR \
 --bqsr-recal-file $path/C113.sorted.markdup.recal_data.table \
 -R $refgenome \
 -I /public1/home/sc30797/ont_data/ngs_bams/aln_bam2/C113.q20.sorted.picard.bam \
 -O $path/C113_marked.sorted.BQSR.bam && echo "** ApplyBQSR done **" &
wait
```

# 5.calling gvcf
GATK4可以允许群体中所有样本都单独与参考基因组call变异，写入.g.vcf.gz文件中取。.g.vcf.gz与.vcf.gz的不同就在于，即便是与参考基因组相同基因型的位点也会被记录下来，以便于整合所有样本。bcftools好像也有这样的功能了，但我没有试过。

perl脚本批量生成calling gvcf命令
```shell
perl creat_hc.pl
```

这一步也可以像bcftools一样将基因组拆成片段来call变异，这回将其拆成10条染色体，这样总命令数就是106 * 10 = 1060，我们依然使用并行提交的方式。

calling gvcf脚本示例
```shell
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa
outdir=~/ont_data/ngs_bams/gatk
$gatk  --java-options "-Xmx4g" HaplotypeCaller \
  --emit-ref-confidence GVCF \
 -R $refgenome \
 -I /public1/home/sc30797/ont_data/ngs_bams/bqsr/C003_marked.sorted.BQSR.bam \
 -L chr1 \
 -O $outdir/C003_marked.HC.chr1.g.vcf.gz && echo "** C003_marked.HC.1.g.vcf.gz done **" &
wait
```

# 6.CombineGVCFs and joint calling
将之前的1060个.g.vcf.gz按照染色体整合起来，形成10个染色体的.g.vcf.gz，再GenotypeGVCFs得到10个染色体的.vcf.gz

perl批量生成gatk命令
```shell
perl create_bash_file_for_Combine.pl
```

gatk脚本示例
```shell
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa
outdir=/public1/home/sc30797/ont_data/ngs_bams/gatk/population
indir=/public1/home/sc30797/ont_data/ngs_bams/gatk
$gatk  --java-options "-Xmx128g" CombineGVCFs  \
 -R $refgenome \
 -V $indir/B73_marked.HC.chr1.g.vcf.gz \
 -V $indir/C002_marked.HC.chr1.g.vcf.gz \
 -V $indir/C003_marked.HC.chr1.g.vcf.gz \
 -V $indir/C004_marked.HC.chr1.g.vcf.gz \
 -V $indir/C005_marked.HC.chr1.g.vcf.gz \
 -V $indir/C006_marked.HC.chr1.g.vcf.gz \
 -V $indir/C007_marked.HC.chr1.g.vcf.gz \
 -V $indir/C008_marked.HC.chr1.g.vcf.gz \
 -V $indir/C009_marked.HC.chr1.g.vcf.gz \
 -V $indir/C010_marked.HC.chr1.g.vcf.gz \
 -V $indir/C011_marked.HC.chr1.g.vcf.gz \
 -V $indir/C012_marked.HC.chr1.g.vcf.gz \
 -V $indir/C013_marked.HC.chr1.g.vcf.gz \
 -V $indir/C014_marked.HC.chr1.g.vcf.gz \
 -V $indir/C017_marked.HC.chr1.g.vcf.gz \
 -V $indir/C018_marked.HC.chr1.g.vcf.gz \
 -V $indir/C019_marked.HC.chr1.g.vcf.gz \
 -V $indir/C020_marked.HC.chr1.g.vcf.gz \
 -V $indir/C021_marked.HC.chr1.g.vcf.gz \
 -V $indir/C022_marked.HC.chr1.g.vcf.gz \
 -V $indir/C023_marked.HC.chr1.g.vcf.gz \
 -V $indir/C024_marked.HC.chr1.g.vcf.gz \
 -V $indir/C025_marked.HC.chr1.g.vcf.gz \
 -V $indir/C026_marked.HC.chr1.g.vcf.gz \
 -V $indir/C028_marked.HC.chr1.g.vcf.gz \
 -V $indir/C029_marked.HC.chr1.g.vcf.gz \
 -V $indir/C032_marked.HC.chr1.g.vcf.gz \
 -V $indir/C033_marked.HC.chr1.g.vcf.gz \
 -V $indir/C034_marked.HC.chr1.g.vcf.gz \
 -V $indir/C035_marked.HC.chr1.g.vcf.gz \
 -V $indir/C036_marked.HC.chr1.g.vcf.gz \
 -V $indir/C037_marked.HC.chr1.g.vcf.gz \
 -V $indir/C038_marked.HC.chr1.g.vcf.gz \
 -V $indir/C039_marked.HC.chr1.g.vcf.gz \
 -V $indir/C041_marked.HC.chr1.g.vcf.gz \
 -V $indir/C043_marked.HC.chr1.g.vcf.gz \
 -V $indir/C044_marked.HC.chr1.g.vcf.gz \
 -V $indir/C046_marked.HC.chr1.g.vcf.gz \
 -V $indir/C047_marked.HC.chr1.g.vcf.gz \
 -V $indir/C048_marked.HC.chr1.g.vcf.gz \
 -V $indir/C049_marked.HC.chr1.g.vcf.gz \
 -V $indir/C050_marked.HC.chr1.g.vcf.gz \
 -V $indir/C051_marked.HC.chr1.g.vcf.gz \
 -V $indir/C052_marked.HC.chr1.g.vcf.gz \
 -V $indir/C053_marked.HC.chr1.g.vcf.gz \
 -V $indir/C054_marked.HC.chr1.g.vcf.gz \
 -V $indir/C057_marked.HC.chr1.g.vcf.gz \
 -V $indir/C058_marked.HC.chr1.g.vcf.gz \
 -V $indir/C059_marked.HC.chr1.g.vcf.gz \
 -V $indir/C060_marked.HC.chr1.g.vcf.gz \
 -V $indir/C066_marked.HC.chr1.g.vcf.gz \
 -V $indir/C067_marked.HC.chr1.g.vcf.gz \
 -V $indir/C068_marked.HC.chr1.g.vcf.gz \
 -V $indir/C069_marked.HC.chr1.g.vcf.gz \
 -V $indir/C070_marked.HC.chr1.g.vcf.gz \
 -V $indir/C073_marked.HC.chr1.g.vcf.gz \
 -V $indir/C074_marked.HC.chr1.g.vcf.gz \
 -V $indir/C075_marked.HC.chr1.g.vcf.gz \
 -V $indir/C078_marked.HC.chr1.g.vcf.gz \
 -V $indir/C080_marked.HC.chr1.g.vcf.gz \
 -V $indir/C081_marked.HC.chr1.g.vcf.gz \
 -V $indir/C082_marked.HC.chr1.g.vcf.gz \
 -V $indir/C083_marked.HC.chr1.g.vcf.gz \
 -V $indir/C084_marked.HC.chr1.g.vcf.gz \
 -V $indir/C086_marked.HC.chr1.g.vcf.gz \
 -V $indir/C087_marked.HC.chr1.g.vcf.gz \
 -V $indir/C088_marked.HC.chr1.g.vcf.gz \
 -V $indir/C089_marked.HC.chr1.g.vcf.gz \
 -V $indir/C090_marked.HC.chr1.g.vcf.gz \
 -V $indir/C092_marked.HC.chr1.g.vcf.gz \
 -V $indir/C094_marked.HC.chr1.g.vcf.gz \
 -V $indir/C095_marked.HC.chr1.g.vcf.gz \
 -V $indir/C096_marked.HC.chr1.g.vcf.gz \
 -V $indir/C097_marked.HC.chr1.g.vcf.gz \
 -V $indir/C099_marked.HC.chr1.g.vcf.gz \
 -V $indir/C100_marked.HC.chr1.g.vcf.gz \
 -V $indir/C101_marked.HC.chr1.g.vcf.gz \
 -V $indir/C103_marked.HC.chr1.g.vcf.gz \
 -V $indir/C104_marked.HC.chr1.g.vcf.gz \
 -V $indir/C105_marked.HC.chr1.g.vcf.gz \
 -V $indir/C106_marked.HC.chr1.g.vcf.gz \
 -V $indir/C107_marked.HC.chr1.g.vcf.gz \
 -V $indir/C108_marked.HC.chr1.g.vcf.gz \
 -V $indir/C109_marked.HC.chr1.g.vcf.gz \
 -V $indir/C110_marked.HC.chr1.g.vcf.gz \
 -V $indir/C111_marked.HC.chr1.g.vcf.gz \
 -V $indir/C112_marked.HC.chr1.g.vcf.gz \
 -V $indir/C113_marked.HC.chr1.g.vcf.gz \
 -V $indir/C114_marked.HC.chr1.g.vcf.gz \
 -V $indir/C115_marked.HC.chr1.g.vcf.gz \
 -V $indir/C121_marked.HC.chr1.g.vcf.gz \
 -V $indir/C122_marked.HC.chr1.g.vcf.gz \
 -V $indir/C124_marked.HC.chr1.g.vcf.gz \
 -V $indir/C126_marked.HC.chr1.g.vcf.gz \
 -V $indir/C128_marked.HC.chr1.g.vcf.gz \
 -V $indir/C129_marked.HC.chr1.g.vcf.gz \
 -V $indir/C130_marked.HC.chr1.g.vcf.gz \
 -V $indir/C131_marked.HC.chr1.g.vcf.gz \
 -V $indir/C153_marked.HC.chr1.g.vcf.gz \
 -V $indir/C154_marked.HC.chr1.g.vcf.gz \
 -V $indir/C155_marked.HC.chr1.g.vcf.gz \
 -V $indir/C156_marked.HC.chr1.g.vcf.gz \
 -V $indir/C157_marked.HC.chr1.g.vcf.gz \
 -V $indir/C159_marked.HC.chr1.g.vcf.gz \
 -V $indir/C161_marked.HC.chr1.g.vcf.gz \
 -V $indir/Mo17_marked.HC.chr1.g.vcf.gz \
 -O $outdir/GATK.106.HC.chr1.g.vcf.gz && echo "** GATK.106.HC.chr1.g.vcf.gz done **" 
$gatk  --java-options "-Xmx128g" GenotypeGVCFs  \
 -R $refgenome \
 -V $outdir/GATK.106.HC.chr1.g.vcf.gz \
 -O $outdir/GATK.106.HC.chr1.vcf.gz && echo "** GATK.106.HC.chr1.vcf.gz done **" 
```

# 7.mergeVCF
就是把10个染色体的.vcf.gz合并起来

```shell
gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk
refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa
outdir=/public1/home/sc30797/ont_data/ngs_bams/gatk/population
indir=/public1/home/sc30797/ont_data/ngs_bams/gatk/population
$gatk MergeVcfs \
 -I $indir/GATK.106.HC.chr1.vcf.gz \
 -I $indir/GATK.106.HC.chr10.vcf.gz \
 -I $indir/GATK.106.HC.chr2.vcf.gz \
 -I $indir/GATK.106.HC.chr3.vcf.gz \
 -I $indir/GATK.106.HC.chr4.vcf.gz \
 -I $indir/GATK.106.HC.chr5.vcf.gz \
 -I $indir/GATK.106.HC.chr6.vcf.gz \
 -I $indir/GATK.106.HC.chr7.vcf.gz \
 -I $indir/GATK.106.HC.chr8.vcf.gz \
 -I $indir/GATK.106.HC.chr9.vcf.gz \
 -O $outdir/GATK.106.HC.vcf.gz && echo "** MergeVcfs done **"
```

# 8.filter vcf
##1.从calling结果中提取到SNP信息

```
/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk SelectVariants -V /public1/home/sc30797/wangzijian/Supplement/GATKout/population2/GATK_113.HC.vcf.gz -select-type SNP -O GATK_113.Wraw-snp.vcf.gz
```

##2.对SNP进行硬过滤
```
/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk VariantFiltration \
    -V GATK_113.Wraw-snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O GATK_113.snps_filtered-2.vcf.gz
```
参数依据：
[Filter variants either with VQSR or by hard-filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

##3.仅保留纯合位点
```
/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk VariantFiltration \
    -V GATK_113.snps_filtered-2.vcf.gz \
    --genotype-filter-expression "isHet == 1" \
    --genotype-filter-name "isHetFilter" \
    -O GATK_113.snps_filtered-hom.vcf.gz

/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk SelectVariants -V GATK_113.snps_filtered-hom.vcf.gz  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -select-type SNP -restrict-alleles-to BIALLELIC -O t113-hard2-snp.vcf.gz
```

##4.提取某个样本的SNP信息
```
/public1/home/sc30797/david/software/gatk-4.1.5.0/gatk SelectVariants  --set-filtered-gt-to-nocall --exclude-filtered --exclude-non-variants  -restrict-alleles-to BIALLELIC  -R Zm-Mo17-REFERENCE-CAU-1.0.fa -V t113-hard2-snp.vcf.gz -sn C78  -O C78-hard2.vcf 
```

