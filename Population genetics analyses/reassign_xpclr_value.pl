#!usr/bin/perl
use strict;
use warnings;

my($file)=@ARGV;

open F,$file or die;

my %hash;
while(<F>){
	chomp;
	my @arr=split;
	my $id="$arr[0]\t$arr[1]\t$arr[2]";
	my $value=$arr[-2];
	if(exists $hash{$id}){
		if($hash{$id}<$value){
			$hash{$id}="$value";
		}
	}
	else{
		$hash{$id}="$value";
	}
}

foreach my $key(keys %hash){
	print $key,"\t","merged\tmerged\t",$hash{$key},"\n";
}

