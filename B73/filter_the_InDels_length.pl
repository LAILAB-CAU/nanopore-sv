#!/usr/bin/env perl
use strict;
use warnings;

die usage() if @ARGV == 0;

my ($vcf,$INS,$DEL) = @ARGV;

open NEW,"$vcf" or die;
open NEW1,">$INS" or die;
open NEW2,">$DEL" or die;

while(<NEW>){
	chomp;
	my @array = split /\s+/;
	if(/^#/){
		print NEW1 "$_\n";
		print NEW2 "$_\n";
	}
	elsif($array[4] eq "<INS>"){
		$array[7] =~ /SVLEN=(\d+);/;
		if($1 >= 50){
			print NEW1 "$_\n";
		}
	}
	elsif($array[4] eq "<DEL>"){
		$array[7] =~ /END=(\d+);/;
		my $length = $1 - $array[1];
		if($length >= 50){
			print NEW2 "$_\n";
		}
	}
}

close NEW;
close NEW1;
close NEW2;


sub usage{
	my $die =<<DIE;
	usage : perl *.pl B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.vcf B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.INS.vcf B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.DEL.vcf
DIE
}
