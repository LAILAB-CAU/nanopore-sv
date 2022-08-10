#!/usr/bin/env perl
use strict;
use warnings;

open NEW,"$ARGV[0]" or die;
#open NEW,"../mo17.EDTA.original.TE.type.bed" or die;  # 新版的EDTA结果跟老版的不一样，需要将gff3处理成bed文件再用
open NEW1,">./$ARGV[0].DNA_DTA.fasta" or die;
open NEW2,">./$ARGV[0].DNA_DTC.fasta" or die;
open NEW3,">./$ARGV[0].DNA_DTH.fasta" or die;
open NEW4,">./$ARGV[0].DNA_DTM.fasta" or die;
open NEW5,">./$ARGV[0].DNA_DTT.fasta" or die;
open NEW6,">./$ARGV[0].DNA_Helitron.fasta" or die;
open NEW7,">./$ARGV[0].MITE_DTA.fasta" or die;
open NEW8,">./$ARGV[0].MITE_DTC.fasta" or die;
open NEW9,">./$ARGV[0].MITE_DTH.fasta" or die;
open NEW10,">./$ARGV[0].MITE_DTM.fasta" or die;
open NEW11,">./$ARGV[0].MITE_DTT.fasta" or die;

while(<NEW>){
	chomp;
	next if (/^#/);
	my @array = split /\s+/;
	#my $start = $array[3] - 1;
	#my $end = $array[4];

	if($array[3] eq "DNA/DTA" or $array[3] eq "DNA/DTC" or $array[3] eq "DNA/DTH" or $array[3] eq "DNA/DTM" or $array[3] eq "DNA/DTT" or $array[3] eq "DNA/Helitron" or $array[3] eq "MITE/DTA" or $array[3] eq "MITE/DTC" or $array[3] eq "MITE/DTH" or $array[3] eq "MITE/DTM" or $array[3] eq "MITE/DTT"){
		my $start = $array[1] - 1;
		my $end = $array[2];

		system("twoBitToFa -seq=$array[0] -start=$start -end=$end Zm-Mo17-REFERENCE-CAU-1.0.2bit ./$ARGV[0].temp.fasta");
		open AA,"./$ARGV[0].temp.fasta" or die;
		while(<AA>){
			chomp;
			if(/^>/){
				if($array[3] eq "DNA/DTA"){
					print NEW1 ">DNA_DTA_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "DNA/DTC"){
					print NEW2 ">DNA_DTC_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "DNA/DTH"){
					print NEW3 ">DNA_DTH_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "DNA/DTM"){
					print NEW4 ">DNA_DTM_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "DNA/DTT"){
					print NEW5 ">DNA_DTT_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "DNA/Helitron"){
					print NEW6 ">DNA_Helitron_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "MITE/DTA"){
					print NEW7 ">MITE_DTA_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "MITE/DTC"){
					print NEW8 ">MITE_DTC_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "MITE/DTH"){
					print NEW9 ">MITE_DTH_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "MITE/DTM"){
					print NEW10 ">MITE_DTM_$array[0]_$start\_$end\n";
				}
				elsif($array[3] eq "MITE/DTT"){
					print NEW11 ">MITE_DTT_$array[0]_$start\_$end\n";
				}
			}
			else{
				if($array[3] eq "DNA/DTA"){
					print NEW1 "$_\n";
				}
				elsif($array[3] eq "DNA/DTC"){
					print NEW2 "$_\n";
				}
				elsif($array[3] eq "DNA/DTH"){
					print NEW3 "$_\n";
				}
				elsif($array[3] eq "DNA/DTM"){
					print NEW4 "$_\n";
				}
				elsif($array[3] eq "DNA/DTT"){
					print NEW5 "$_\n";
				}
				elsif($array[3] eq "DNA/Helitron"){
					print NEW6 "$_\n";
				}
				elsif($array[3] eq "MITE/DTA"){
					print NEW7 "$_\n";
				}
				elsif($array[3] eq "MITE/DTC"){
					print NEW8 "$_\n";
				}
				elsif($array[3] eq "MITE/DTH"){
					print NEW9 "$_\n";
				}
				elsif($array[3] eq "MITE/DTM"){
					print NEW10 "$_\n";
				}
				elsif($array[3] eq "MITE/DTT"){
					print NEW11 "$_\n";
				}
			}
		}
		close AA;
		system("rm ./$ARGV[0].temp.fasta");
	}
}

close NEW;
close NEW1;
close NEW2;
close NEW3;
close NEW4;
close NEW5;
close NEW6;
close NEW7;
close NEW8;
close NEW9;
close NEW10;
close NEW11;
