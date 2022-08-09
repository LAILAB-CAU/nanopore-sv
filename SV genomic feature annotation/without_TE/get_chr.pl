#!/usr/bin/perl
use strict;
use warnings;

my($file,$chr)=@ARGV;

open F,$file or  die;
while(<F>){
	chomp;
	my @arr=split;
	if($arr[0] eq "chr$chr"){
		print $_,"\n";
	}
}
