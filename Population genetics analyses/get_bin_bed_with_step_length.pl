#!/usr/bin/perl
use strict;
use warnings;

my($file,$size,$step)=@ARGV;
open F,$file or die;

my %hash;
my %num;
while(<F>){
	chomp;
	my @arr=split;
	my $chr=$arr[0];
	my $length=$arr[1];
	$hash{$chr}=$length;
}

for(my $i=1;$i<=10;$i++){
	my $chr="chr$i";
	my $length=$hash{$chr};
	for(my $j=1;$j<=$length;$j+=$step){
		my $start=$j;
		my $end=$start+$size-1;
		if($end<$length){
			print $chr,"\t",$start,"\t",$end,"\n";
		}
		else{
			$end=$length;
			print $chr,"\t",$start,"\t",$end,"\n";
			last;
		}
	}
}
			

