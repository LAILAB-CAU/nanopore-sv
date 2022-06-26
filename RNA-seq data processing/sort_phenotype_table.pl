#!/usr/bin/perl
use strict;
use warnings;

my($file,$ref1,$ref2)=@ARGV;

open F,$file or die;
open R1,$ref1 or die;
open R2,$ref2 or die;

my $i=0;
my %hash;
while(<R2>){
	chomp;
	$i++;
	my @arr=split;
	my $id=$arr[1];
	$hash{$i}=$id;
}

my %name;
while(<R1>){
	chomp;
	my @arr=split;
	for(my $j=1;$j<@arr;$j++){
		$name{$arr[$j]}=$j;
	}
	last;
}

while(<F>){
	chomp;
	my @arr=split;
	my $seq=$arr[0];
	for(my $i=1;$i<@arr;$i++){
		my $info=$name{$hash{$i}};
		$seq.="\t$arr[$info]";
	}
	print $seq,"\n";
}

	
