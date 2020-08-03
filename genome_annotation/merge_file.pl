#!/usr/bin/perl
use strict;

my($file)=@ARGV;
open F,$file or die;

my $start=-100;
my $end=-100;
my $chr;
while(<F>){
	chomp;
	my @arr=split;
	if($arr[1]==($end+1) && $chr eq "$arr[0]"){
		$end=$arr[1];
	}
	else{
		unless($start==-100){
			
			print $chr,"\t",$start,"\t",$end,"\n";
		}
		$start=$arr[1];
		$end=$arr[1];
		$chr=$arr[0];
	}
}

unless($start==-100){
	print $chr,"\t",$start,"\t",$end,"\n";
                }

