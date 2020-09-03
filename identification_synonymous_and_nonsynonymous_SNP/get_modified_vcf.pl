#!/usr/bin/perl
use strict;

my($file,$file_out)=@ARGV;

open F,$file or die;
open O,">",$file_out or die;

while(<F>){
	chomp;
	my @arr=split;
	if(/#/){
		next;
	}
	print O $arr[0],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t",$arr[4],"\n";
}


