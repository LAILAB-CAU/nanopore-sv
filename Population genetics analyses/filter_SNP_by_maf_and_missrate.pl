#!/usr/bin/perl
use strict;
use warnings;

my($file,$file_out)=@ARGV;

open F,"zcat $file|" or die;
open O,"|bgzip  > $file_out" or die;
while(<F>){
	chomp;
	my @arr = split /\s+/;
	if(/^#/){
		print O $_,"\n";
	}
	else{
		my $missing_number = 0;
		my $ref_number = 0;
		my $alt_number = 0;
		my $total;
		for(my $number = 9;$number < @arr;$number++){
			my @temp=split /:/,$arr[$number];
			$arr[$number]=$temp[0];
			if($arr[$number] eq "0/1" or $arr[$number] eq "1/0"){
				$ref_number++;
				$alt_number++;
				$total++;
			}
			elsif( $arr[$number] eq "./." ){
				$missing_number++;
			}
			elsif($arr[$number] eq "0/0"){
				$ref_number+=2;
				$total++;
			}
			elsif($arr[$number] eq "1/1"){
				$alt_number+=2;
				$total++;
			}
		}

		my $missing = $missing_number/$total;
		next if ($missing >= 0.5);
		my $allele=$ref_number + $alt_number;
		if($ref_number/$allele >=0.05 && $alt_number/$allele >=0.05){
			print O "$_\n";
		}
	}
}
close F;

