#!/usr/bin/env perl
use strict;
use warnings;

my ($file_path) = @ARGV;
my %hash_sample_fq1;
my %hash_sample_fq2;
open NEW, "$file_path" or die;
while(<NEW>){
    chomp;
    my @a = split /\s+/;
    $hash_sample_fq1{$a[0]} = $a[1];
    $hash_sample_fq2{$a[0]} = $a[2];
}
close NEW;


open NEW2, ">sampe" or die;
foreach my $sample(keys %hash_sample_fq1){
    open NEW,">./$sample\_bwa_q20_sort.sh" or die;
    print NEW "#!/bin/bash\n";
    print NEW "refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa\n";

    #for already fastp fq
    # print NEW "#bwa aln -t 32  \$refgenome  $hash_sample_fq1{$sample} -f $sample\_1.sai &\n";
    # print NEW "#bwa aln -t 32  \$refgenome  $hash_sample_fq2{$sample} -f $sample\_2.sai &\n";
    # print NEW "#wait\n";

    #for raw fq
    print NEW "fastp -i $hash_sample_fq1{$sample}  -I $hash_sample_fq2{$sample}  -o $sample\_fastp1.fq.gz  -O $sample\_fastp2.fq.gz  -j $sample.json  -h $sample.html -w 16\n";
    print NEW "bwa aln -t 32  \$refgenome  $sample\_fastp1.fq.gz -f $sample\_1.sai &\n";
    print NEW "bwa aln -t 32  \$refgenome  $sample\_fastp2.fq.gz -f $sample\_2.sai &\n";
    print NEW "wait\n";

    print NEW2 "bwa sampe -r \"\@RG\\tID:$sample\\tPL:ILLUMINA\\tLB:$sample\\tSM:$sample\"  ";
    #for already fastp fq
    #print NEW2 "\$refgenome $sample\_1.sai $sample\_2.sai $hash_sample_fq1{$sample} $hash_sample_fq2{$sample}  ";
    #for raw fq
    print NEW2 "\$refgenome $sample\_1.sai $sample\_2.sai $sample\_fastp1.fq.gz $sample\_fastp2.fq.gz  ";
    print NEW2 "\| samtools view -b   -q 20 -  ";
    print NEW2 "\| samtools sort   - > $sample.q20.sorted.bam\n";
    print NEW2 "samtools index $sample.q20.sorted.bam\n ";
    close NEW;
}
close NEW2;
