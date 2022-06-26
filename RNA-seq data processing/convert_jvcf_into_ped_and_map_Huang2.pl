#!/usr/bin/perl

use strict;
use warnings;
die usage() if @ARGV == 0;
my ($jvcf,$ped,$map) = @ARGV;

open P,">$ped" or die;
open M,">$map" or die;
open F,"zcat $jvcf | " or die;

my %hash_genotype;
my %hash_accessions;
while(<F>){
	chomp;
	my @array = split /\s+/;
	if(/^#CHROM/){
		for(my $number = 3;$number < @array;$number++){
			$hash_accessions{$number} = $array[$number];
		}
	}
	else{
		$array[0] =~ s/chr//g;
		print M "$array[0]\t$array[0]_$array[1]\t$array[1]\n";
		
		my $ref;
		my $alt;
		if($array[2]=~/(.)\/(.)/){
			$ref=$1;
			$alt=$2;
		}
		else{
			die;
		}
		for(my $number = 3;$number < @array;$number++){
			if($array[$number] eq "X" ){
				$hash_genotype{$number}.="\t0 0";
			}
			elsif($array[$number] eq "x"){
				$hash_genotype{$number}.="\t$ref $alt";
			}
			else{
				$hash_genotype{$number}.="\t$array[$number] $array[$number]";
			}
		}
	}

}

close F;

my $family_cal;
foreach my $key(sort {$a <=> $b} keys %hash_accessions){
	$family_cal++;
	print P "FAM$family_cal\t$hash_accessions{$key}\t0\t0\t0\t0$hash_genotype{$key}\n";
}

close P;
close M;

sub usage{
	my $die =<<DIE;
	usage : perl *.pl *.jvcf *.ped *.map
DIE
}
	
