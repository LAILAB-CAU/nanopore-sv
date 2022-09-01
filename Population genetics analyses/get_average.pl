#!/usr/bin/perl
use strict;
use warnings;

my($file,$ref)=@ARGV;

open F,$file or die;
open R,$ref or die;

my %bin;
my %hash;
while(<F>){
	chomp;
	my @arr=split;
	my $id="$arr[0]\t$arr[1]\t$arr[2]";
	if( $arr[4]==-1){
		next;
	}
	else{
		if($arr[6]<10){
			next;
		}
		$hash{$id}+=$arr[-2];
		$bin{$id}++;
	}
}

while(<R>){
	chomp;
	my @arr=split;
	my $id="$arr[0]\t$arr[1]\t$arr[2]";
	my $ave;	
	if(exists $hash{$id}){
		$ave=$hash{$id}/$bin{$id};
	}
	else{
		$hash{$id}="NA";
		$bin{$id}="0";
		$ave="NA";
	}
	print $id,"\t",$hash{$id},"\t",$bin{$id},"\t",$ave,"\n";
}
	
