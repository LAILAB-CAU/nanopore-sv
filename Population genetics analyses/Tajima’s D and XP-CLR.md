Tajima's D was calculated for each 10 kbp genomic window by the vcftools (v 0.1.16, --TajimaD 10000), based on SNPs with missing rate less than 0.5. 

```
DEL=DEL.SViper.out.1000.vcf.uniqe.nomissing.vcf
INS=INS.SViper.out.1000.vcf.uniqe.nomissing.vcf
vcftools --vcf $INS --TajimaD 10000 --out 106_INS_TajiamD
vcftools --vcf $DEL --TajimaD 10000 --out 106_DEL_TajiamD


SNP=GATK.106.snp.after_hard_filtered.vcf.gz
vcftools --gzvcf $SNP --TajimaD 10000 --out 106_SNP_TajiamD2 --max-missing 0.5 
```

The XP-CLR values for 10 kbp windows were calculated by XP-CLR (https://github.com/hardingnj/xpclr) (--ld 0.95 --maxsnps 200 --size 10000 --step 1000), based on with missing rate less than 0.5 and MAF larger than 0.05. The XP-CLR values were calculated between maize groups with top and bottom 20% PH. Adjacent 10 kbp windows from the top 20% of the XP-CLR results were merged into larger regions and only one window lower than 20% was retained for each region. Regions in the highest 5% of the mean region-wise XP-CLR values were regarded as selective regions. The genome-wide cutoff for XP-CLR was set to as the lower bound of the top 5% XP-CLR values.

```
perl filter_SNP_by_maf_and_missrate.pl chr1.raw.vcf2.gz chr1.vcf.gz 
tabix chr1.vcf.gz
xpclr --format vcf --samplesA PH_H.list2 --samplesB PH_L.list2 --chr chr1 --input $fold/chr1.vcf.gz   --ld 0.95 --maxsnps 200 --size 10000 --step 1000  --out ./PH_chr1.txt &
perl get_NA_zero_xpclr_bin_num.pl ../PH_chr1.txt chr/chr1_xpclr.bed

perl get_bin_bed_with_step_length.pl chrom_length 10000 10000 > B73_10kb_bin_step10k
cat chr/chr*_xpclr.bed > chr/whole_genome_xpclr.bed

bedtools intersect -a B73_10kb_bin_step10k -b chr/whole_genome_xpclr.bed  -wao  |perl get_average.pl  - B73_10kb_bin_step10k > whole_genome_averaged_10kb_xpclr

perl get_the_significant_windows_v2.pl whole_genome_averaged_10kb_xpclr 0.2 whole_genome_averaged_10kb_xpclr_top20
bedtools merge -i whole_genome_averaged_10kb_xpclr_top20 -d 10000 |bedtools intersect -a - -b whole_genome_averaged_10kb_xpclr_top20 -wao  |perl reassign_xpclr_value.pl  -  > whole_genome_averaged_10kb_xpclr_top20_merge
bedtools  intersect -a whole_genome_averaged_10kb_xpclr -b whole_genome_averaged_10kb_xpclr_top20_merge -v |cat - whole_genome_averaged_10kb_xpclr_top20_merge > whole_genome_region_wise_xpclr
perl get_the_significant_windows_v2.pl whole_genome_region_wise_xpclr 0.05 whole_genome_region_wise_xpclr_top_0.05
```