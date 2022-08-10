#!/usr/bin/perl
use strict;
my($file,$file_out,$id_out,$k)=@ARGV;
open F,$file or die;
open O,">",$file_out or die;
open ID,">",$id_out or die;
my %num;
my $num;
my %length;
my %chr_name;
my %pos;
my $id=0;
my %seq;
while(<F>){
	chomp;
	if(/##/){
		next;
	}
	if(/CHROM/)	{
		next;
	}
	my @arr=split;
	$id++;
	my $length=-1;
#	length=length($arr[3])-length($arr[4]);
	my @temp=split /;/,$arr[7];
	if($temp[$k]=~ s/SVLEN=//g){
		if($temp[$k] <0){
			$temp[$k]=-$temp[$k];
		}
		$length=$temp[$k];
	}
	else{
		die;
	}
	#
	if($length<0){
		$length=0-$length;
	}	
	$num{$id}++;
	$length{$id}=$length;
	$chr_name{$id}=$arr[0];
	$pos{$id}=$arr[1];
	$seq{$id}=$arr[4];
}
my $total;
foreach my $key (keys %num){
	$total++;
	print O ">",$key,"\n";
	print O  $seq{$key},"\n";
	print ID $key,"\t",$chr_name{$key},"_", $pos{$key},"_",$pos{$key}+$length{$key}, "\n";
}
print $file,"\t",$total,"\n";