# 1. Base calling

The raw file format produced by the ONT PromethION sequencer is a binary fast5 format. Runs were base-called with Guppy (v 3.3.0) (https://community.nanoporetech.com) with the PromethION r9.4.1 high accuracy (HAC) model.
```
#!/bin/bash
ont-guppy/bin/guppy_basecaller --version           >>./version-check.log
ont-guppy/bin/guppy_basecaller  --print_workflows  >>./config-check.log
ont-guppy/bin/guppy_basecaller	--device cuda:0,1 --num_callers 10   --ipc_threads 20 --gpu_runners_per_device 4  --chunks_per_runner 1664	--flowcell	FLO-PRO002	--kit	SQK-LSK109	--recursive	--qscore_filtering	--verbose_logs	--input_path	/fast5_path	--save_path	/xxx
```
# 2. Self-correction, read mapping and preliminary SV calling

Take B073 as an example: 
```
The raw sequencing fastq file name (After base calling): B073.hac.raw.fastq
Mo17 genome reference file name: Maize_Mo17.fasta (downloaded from https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/)
```
Step 1. Generate PAF file to obtain the statistics of the raw data mapping by using minimap2 (https://github.com/lh3/minimap2):
```
minimap2 -c -M 0 --secondary=no -Q -t 10 -x map-ont Maize_Mo17.fasta B073.hac.raw.fastq > B073toMo17.raw.paf
sh statistics_paf.sh B073toMo17.raw.paf
```
Step 2. Raw data correction using NECAT version 1 (https://github.com/xiaochuanle/NECAT) by using the Correct module with default parameters:
```  
perl necat.pl correct run.B073.cfg
```
The name of corrected reads file: B073.hac.cns.fasta.

Step 3. Generate PAF file to obtain the statistics of the corrected data mapping: 
```
minimap2 -c -M 0 --secondary=no -Q -t 10 -x map-ont Maize_Mo17.fasta B073.hac.cns.fastq > B073toMo17.cns.paf
sh statistics_paf.sh B073toMo17.cns.paf
```
Step 4. SV calling using PBSV (https://github.com/PacificBiosciences/pbsv):
```
pbmm2 align --sort -j 31 -J 8 -m 1000M --rg '@RG\tID:B073\tSM:B073' Maize_Mo17.fasta B073.hac.cns.fasta B073.bam
pbsv discover -q 10 B073.bam B073.svsig.gz
pbsv call -j 10 -m 20 --min-cnv-length 1K --max-ins-length 10K --max-dup-length 100K Maize_Mo17.fasta B073.svsig.gz B073.PBSV.vcf
```
Step 5. Filte PBSV results by applying the following crieteria:

(1) Only consider Deletion and Insertion 

(2) Only consider the genotype: 1/1

(3) The depth of reference allele: 0 

(4) Total depth for a SV: at least 2
```
perl sel_pbsv_raw_format_vcf.pl B073.PBSV.vcf B073.PBSV.vcf.filter
```
Step 6. Using the Bam file from pbmm2 to recall SV using Sniffles software (https://github.com/fritzsedlazeck/Sniffles): 
```
# add MD field: 
  samtools calmd -@ 8 B073.bam Maize_Mo17.fasta | samtools sort - -@ 8 -m 1G -o B073.MD.bam
# Run SNIFFLES (V1.0.11): 
 sniffles -m B073.MD.bam -v B073.sniff.vcf -s 2 -t 4 -l 20 --report_BND --report_seq --report_read_strands --ignore_sd --tmp_file ./tmpfile/B073.tmp --cluster --cluster_support 1 --genotype
```
Step 7. Filter Sniffles results (using the same criterion with the results of PBSV):
```
perl sel_sniffles_raw_format_vcf.pl B073.sniff.vcf B073.sniff.vcf.filter
```

Step 8. For each accession the final set is the intersect of both PBSV and Sniffles results. This was done by SURVIVOR (https://github.com/fritzsedlazeck/SURVIVOR): 
```
SURVIVOR merge B073_DEL_PBSV_SNIFFLES_filtered.lst 1000 1 1 -1 -1 -1 B073.DEL.intersect.1000.vcf
SURVIVOR merge B073_INS_PBSV_SNIFFLES_filtered.lst 1000 1 1 -1 -1 -1 B073.INS.intersect.1000.vcf
```

# 3. SV breakpoint refinement with SRS data

Step 1. Prepare short-read sequencing mapping bam files for each accession, see the "Short-read sequencing data processing" section. Suppose the bam file for B073 is B073.sorted.bam.

Step 2. SViper installation:

SViper (https://github.com/smehringer/SViper) polishes deletions and insertions called on long read data (ONT) using short exact reads for refinement. Using --recursive mode when downloading the package through git clone so that the third-party extern package could be downlowded as well, otherwise, errors will occur during installation.
  
Step 3. Running SViper to polish deletions and insertions:

B073.bam is the mapping file for ONT long-read data. 

```
$SViper_folder/build/sviper -c B073.DEL.intersect.1000.vcf -s B073.sorted.bam -l B073.bam -r Zm-Mo17-REFERENCE-CAU-1.0.fa -t 60 -o B073.SViper.DEL
$SViper_folder/build/sviper -c B073.INS.intersect.1000.vcf -s B073.sorted.bam -l B073.bam -r Zm-Mo17-REFERENCE-CAU-1.0.fa -t 60 -o B073.SViper.INS
```
     
# 4. Obtaining population SV sets

Combine the polished deletion and insertion VCF into a merged VCF file using SURVIVOR and obtain the population SV sets. Put the file names of polished SV set for all accessions in a list: DEL.combine.1000.lst and INS.combine.1000.lst. 

```
SURVIVOR merge DEL.combine.1000.lst 1000 1 1 -1 -1 -1 DEL.SViper.out.1000.vcf
SURVIVOR merge INS.combine.1000.lst 1000 1 1 -1 -1 -1 INS.SViper.out.1000.vcf
```

# 5. Genotyping. Seperate 0/0 from ./. in the vcf, and filter SVs with high missing rate


Genotyping with ONT bam file. Create a depth file for each ONT data, put all *.depth.gz file in a folder, which will be the input for genotyping.

```
cat DEL.SViper.out.1000.vcf INS.SViper.out.1000.vcf | sortBed - > INDEL.SViper.out.1000.vcf.sorted.bed
srun -n 1 samtools depth -b INDEL.SViper.out.1000.vcf.sorted.bed Mo17.bam | gzip -c - > Mo17.depth.gz

perl genotyping_and_unique.pl DEL.SViper.out.1000.vcf INS.SViper.out.1000.vcf DEL.SViper.out.1000.vcf.uniqe.nomissing.vcf INS.SViper.out.1000.vcf.uniqe.nomissing.vcf
```

The population SV set for deletions and insertions (length > 50bp and unique) are: DEL.SViper.out.1000.vcf.uniqe.nomissing.vcf INS.SViper.out.1000.vcf.uniqe.nomissing.vcf. 








