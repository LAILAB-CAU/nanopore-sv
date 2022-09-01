#!/usr/bin/perl
use strict;
use warnings;

my($file,$file_out)=@ARGV;

open F,$file or die;
open O,">",$file_out or die;

my $num1;
my $num2;
my $num3;

while(<F>){
	chomp;
	if(/xpclr/){
		next;
	}
	$num1++;
	my @arr=split /\t+/;
	unless(@arr>10){
		next;
	}
	$num2++;
	print O "$arr[1]\t$arr[2]\t$arr[3]\t$arr[-3]\t$arr[-2]\n";
	unless($arr[11]>0){
		next;
	}
	$num3++;
}

print "$file\t$num1\t$num2\t$num3\n"	;

