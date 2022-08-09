#!/usr/bin/perl
use strict;
use warnings;

my($file,$file_out)=@ARGV;
open F ,$file or die;
open O ,">",$file_out or die;
while(<F>){
	chomp;
	my @arr=split;
	unless($arr[2] =~ /UTR/ && $arr[2] =~ /five/){
		next;
	}
		my @tempArr=split /;/ ,$arr[8];
	my $parent;
	if($tempArr[0] =~ /Parent=(.*)/){
		$parent=$1;
	}
	else{
		print $_,"\n";
	}
	print O $arr[0],"\t",$arr[3],"\t",$arr[4],"\t",$parent,"\t",$arr[6],"\n";
}

