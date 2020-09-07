#!/usr/bin/env perl
use strict;
use warnings;
my @files = </public1/home/sc30797/ont_data/ngs_bams/bqsr/*_marked.sorted.BQSR.bam>;


my $cal;
my $file_cal;
my @samples;

foreach my $file(@files){
	my @a = split /\//,$file;
	my @b = split /_/,$a[-1];	
	my $sample = $b[0];
    push @samples,$b[0];
}

for(my $i = 1;$i <=10;$i++){
    open NEW,">./chr$i\_Combine_Genotype.sh" or die;
    print NEW "#!/bin/bash\n";
    print NEW "gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk\n";
    print NEW "refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa\n";
    print NEW "outdir=/public1/home/sc30797/ont_data/ngs_bams/gatk/population\n";
    print NEW "indir=/public1/home/sc30797/ont_data/ngs_bams/gatk\n";
    print NEW "\$gatk  --java-options \"-Xmx128g\" CombineGVCFs  \\\n -R \$refgenome \\\n";
    foreach my $very_sample(@samples){
        print NEW " -V \$indir/$very_sample\_marked.HC.chr$i.g.vcf.gz \\\n";
    }
    print NEW " -O \$outdir/GATK.106.HC.chr$i.g.vcf.gz && echo \"** GATK.106.HC.chr$i.g.vcf.gz done **\" \n";

    print NEW "\$gatk  --java-options \"-Xmx128g\" GenotypeGVCFs  \\\n -R \$refgenome \\\n -V \$outdir/GATK.106.HC.chr$i.g.vcf.gz \\\n -O \$outdir/GATK.106.HC.chr$i.vcf.gz && echo \"** GATK.106.HC.chr$i.vcf.gz done **\" \n";
}


