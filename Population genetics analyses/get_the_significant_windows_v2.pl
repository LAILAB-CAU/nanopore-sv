#!/usr/bin/perl
use strict;
use warnings;

my ($file,$significant_cutoff,$output_regions) = @ARGV;

open NEW ,$file or die;
my @array_xpclr;

while(<NEW>){
	chomp;
	my @array = split /\s+/;
	if(/NA/){
		next;
	}
	push @array_xpclr,$array[-1];
}
########## sorted the XPCLR from highest to lowest
my @sorted_array_xpclr = sort {$b <=> $a} @array_xpclr;
my $a = @sorted_array_xpclr;
my $cut = int($a*$significant_cutoff);
my $xpclr_cut = $sorted_array_xpclr[$cut];
print "##################################\nThe XPCLR cutoff under the significance ($significant_cutoff) is $xpclr_cut\n";
open NEW1,">$output_regions" or die;
open NEW ,$file or die;
while(<NEW>){
	chomp;
	my @array = split /\s+/;
	if(/NA/){
		next;
	}
	if($array[-1] >= $xpclr_cut){
		print NEW1 "$_\n";
	}
}
close NEW;
close NEW1;


