#!/usr/bin/perl
use strict;
use warnings;

my($file,$file_out)=@ARGV;
open F ,$file or die;
open O ,">",$file_out or die;

while(<F>){
	chomp;
	my @arr=split;
	unless($arr[2] =~ /gene/){
		next;
	}
	my @tempArr=split /;/ ,$arr[8];
	my $ID;
	if($tempArr[0] =~ /ID=(.*)/){
		$ID=$1;
	}
		else{
		print $_,"\n";
	}
	print O $arr[0],"\t",$arr[3],"\t",$arr[4],"\t",$ID,"\t",$arr[6],"\n";
}


