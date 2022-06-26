#!/usr/bin/perl
use strict;
use warnings;


die usage() if @ARGV == 0;
my ($vcf,$jvcf) = @ARGV;

open F,"zcat $vcf|" or die;
open O,"|gzip  > $jvcf"  or die;

my %hash;
while(<F>){
	chomp;
	next if $_=~ /^##/;
	my @col = split /\s+/;
	if($col[0] eq "#CHROM"){
		next;
	}
	else{
		my $id="$col[0]-$col[1]";
		$hash{$id}++;
	}
}

open F,"zcat $vcf|" or die;
while(<F>){
	chomp;
	next if $_=~ /^##/;
	my @col = split /\s+/;
	if($col[0] eq "#CHROM"){
		print O "#CHROM\tPOS\tAllele\t";
		for(my $col_number = 9;$col_number < @col;$col_number++){
			print O "$col[$col_number]\t";
		}
		print O "\n";
	}
	else{
		my $id="$col[0]-$col[1]";
		if($hash{$id}>1){
			next;
		}

		my $ref_len=length($col[3]);
		my $alt_len=length($col[4]);
		$col[4]="G";
		$col[3]="A";

		print O "$col[0]\t$col[1]\t$col[3]\/$col[4]\t";
		for(my $col_number = 9;$col_number < @col;$col_number++){
			my @snp_infor = split /:/,$col[$col_number];
			my $genotype= $snp_infor[0];
			if($genotype eq '0/0'|| $genotype eq '0|0'){
				 print O "$col[3]\t";
			}
			elsif($genotype eq '0/1'|| $genotype eq '0|1'){
				print O "x\t";
			}
			elsif($genotype eq '1/0'|| $genotype eq '1|0'){
				print O "x\t";
			}
			elsif($genotype eq '1/1'|| $genotype eq '1|1'){
				print O "$col[4]\t";
			}
			elsif($genotype eq './.'|| $genotype eq '.|.' || $genotype eq '.'){
				print O "X\t";
			}
			else{
				print "Error:$col[0]\t$col[1]\t$col[$col_number]\n";
			}
		}
		print O "\n";
	}

}

