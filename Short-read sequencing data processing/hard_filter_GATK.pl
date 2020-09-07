#!/usr/bin/env perl
use strict;
use warnings;

my ($input,$output_vcffile) = @ARGV;
#打开region list文件
#open VCF, "gzip -dc $input |" or die;
open VCF, "$input" or die;

open OUTVCF, ">$output_vcffile" or die;


while(my $bcfline = <VCF>){
    chomp $bcfline;
    if($bcfline =~ /^#/){
        print OUTVCF "$bcfline\n";
        next;
    }
    #读取第10列以后每一列(样本基因型)的信息
    my @vcf_array = split /\s+/, $bcfline;
    #有两个以上的alt genotype的不要
    next if length($vcf_array[4])>1;
    my $sample_sum = $#vcf_array-9;
    #统计该位点杂合数 0/0数 1/1数 missing数
    my $heter_sum = () = $bcfline =~ /0\/1|1\/0/gi;
    my $ref_sum = () = $bcfline =~ /0\/0/gi;
    my $variant_sum = () = $bcfline =~ /1\/1/gi;
    my $missing_sum = () = $bcfline =~ /\.\/\./gi;
    # for(my $i = 9;$i <= $#vcf_array;$i++){
    #     my @genotype = split /:/,$vcf_array[$i];
    #     if($genotype[0] eq '0/0'){$ref_sum++;}
    #     elsif($genotype[0] eq '1/1'){$variant_sum++;}
    #     elsif($genotype[0] eq './.'){$missing_sum++;}
    #     else{$heter_sum++;}
    # }
    #计算missing rate maf 杂合率 
    my $missing_rate = $missing_sum/$sample_sum;
    my $maf;
    next if ($ref_sum+$variant_sum == 0);
    if($ref_sum >= $variant_sum){$maf = $variant_sum/($ref_sum+$variant_sum);}
    else{$maf = $ref_sum/($ref_sum+$variant_sum);}
    my $hetero_rate = $heter_sum/$sample_sum;
    #只输出missing rate低于0.5 maf高于0.05 杂合率低于0.05的位点
    if(($missing_rate<=0.5) and ($maf > 0.05) and ($hetero_rate <= 0.05)){
        print OUTVCF "$bcfline\n";
    }
}
close VCF;
close OUTVCF;

