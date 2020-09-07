#!/usr/bin/env perl
use strict;
use warnings;

die if @ARGV == 0;
open NEW, "$ARGV[0]" or die;
open OUT, ">$ARGV[1]" or die;
my $ref_genome = "/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa";
my $bam_list = "./bamlist";

my $heading = "#!/bin/bash\nbcftools=~/wangzijian/software/bcftools-1.10.2/my_bcftools/bin/bcftools\nexport LD_LIBRARY_PATH=/public1/home/sc30797/ont_data/ngs_bams/bcf_snp/lib:\$LD_LIBRARY_PATH\n";
my $commamd1 = "\$bcftools mpileup -E -C50 -Q20 -q20  -Ou  --threads 21 -r ";
my $command2= " -f $ref_genome -b $bam_list | bcftools call -cv -Oz --threads 21 -o  bcftools_106_";

my $count=0;
while(<NEW>){
	chomp;
    $count+=1;
    my $re = $_ =~ s/[:-]/_/r;
    if($count % 3 == 1){
        print OUT "$heading";
        print OUT "$commamd1 $_ $command2\_$re.vcf.gz &\n";
    }
    else{
        print OUT "$commamd1 $_ $command2\_$re.vcf.gz &\n";
    }
    if($count % 3 == 0 or $count==46){
        print OUT "wait\n";
    }
    
	
}


close NEW;
close OUT;
