#!/usr/bin/perl
use strict;
my($file,$file_out1,$file_out2)=@ARGV;
open F,$file or die;
open O1,">",$file_out1 or die;
open O2,">",$file_out2 or die;

while(<F>){
	chomp;
	if(/##/){
		next;
	}
	my $line=$_;
	$line =~ s/ANN=//g;
	my @TEMP=split /\s+/,$line;
	if($TEMP[4]=~/,/){
		next;
	}
	my @T=split /,/,$TEMP[7];
	#print $TEMP[7];
	if(/missense_variant/){
		print O1 $_,"\n";
	}
	elsif(/synonymous_variant/){
		print O2 $_,"\n";
        }
}
