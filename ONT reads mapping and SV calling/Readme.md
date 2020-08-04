Take B073 as an example: 
```
The raw sequencing fastq file name (Hac mode): B073.hac.raw.fastq
Mo17 genome reference file name:               Maize_Mo17.fasta
```
(1) Generate PAF file to statistics the raw data mapping information: 
```
minimap2 -c -M 0 --secondary=no -Q -t 10 -x map-ont Maize_Mo17.fasta B073.hac.raw.fastq > B073toMo17.raw.paf
sh statistics_paf.sh B073toMo17.raw.paf
```
(2) Raw data correction using Necat version 1:
```  
  perl necat.pl correct run.B073.cfg
  The name of corrected reads file: B073.hac.cns.fasta
```
(3) Generate PAF file to statistics the corrected data mapping information: 
```
minimap2 -c -M 0 --secondary=no -Q -t 10 -x map-ont Maize_Mo17.fasta B073.hac.cns.fastq > B073toMo17.cns.paf
sh statistics_paf.sh B073toMo17.cns.paf
```
(4) SV calling using PBSV
```
1). pbmm2 align --sort -j 31 -J 8 -m 1000M --rg '@RG\tID:B073\tSM:B073' Maize_Mo17.fasta B073.hac.cns.fasta B073.bam
2). pbsv discover -q 10 B073.bam B073.svsig.gz
3). pbsv call -j 10 -m 20 --min-cnv-length 1K --max-ins-length 10K --max-dup-length 100K Maize_Mo17.fasta B073.svsig.gz B073.PBSV.vcf
```
(5). Filte PBSV results, 
```
criterion: 
(1) Only consider Deletion and Insertion 
(2) Only consider the genotype: 1/1
(3) The depth of reference allele: 0 
(4) Total depth for a SV: at least 2
perl sel_pbsv_raw_format_vcf.pl B073.PBSV.vcf B073.PBSV.vcf.filter
```
(6) Using the Bam file from pbmm2 to recall SV using SNIFFLES software: 
```
1) add MD field: 
   samtools calmd -@ 8 B073.bam Maize_Mo17.fasta | samtools sort - -@ 8 -m 1G -o B073.MD.bam
2) Run SNIFFLES (V1.0.11): 
   sniffles -m B073.MD.bam -v B073.sniff.vcf -s 2 -t 4 -l 20 --report_BND --report_seq --report_read_strands --ignore_sd --tmp_file ./tmpfile/B073.tmp --cluster --cluster_support 1 --genotype
```
(7) Filter SNIFFLES results (using the same criterion with the results of PBSV):
```
perl sel_sniffles_raw_format_vcf.pl B073.sniff.vcf B073.sniff.vcf.filter
```
(8) For each sample the final set is the intersect of both PBSV and SNIFFLES results: 
```
SURVIVOR merge B073_DEL_PBSV_SNIFFLES_filtered.lst 1000 1 1 -1 -1 -1 B073_DEL_intersect.vcf
SURVIVOR merge B073_INS_PBSV_SNIFFLES_filtered.lst 1000 1 1 -1 -1 -1 B073_INS_intersect.vcf
```
(9) Using the intersect result of each sample, combine them in a merged VCF file using SURVIVOR:
```
SURVIVOR merge DEL_pbsv_sniffles_intersect.lst 1000 1 1 -1 -1 -1 DEL.merge.vcf
SURVIVOR merge INS_pbsv_sniffles_intersect.lst 1000 1 1 -1 -1 -1 INS.merge.vcf
```
