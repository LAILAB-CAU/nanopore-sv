#!/usr/bin/env perl
use strict;
use warnings;

die usage() if @ARGV == 0;

my %hash_sequence;
my $header;
open NEW,"$ARGV[0]" or die;
while(<NEW>){
	chomp;
	if(/>/){
		$_ =~ s/>//;
		$header = $_;
	}
	else{
		$hash_sequence{$header} .= $_;
	}
}
close NEW;

my %hash_cal;
open NEW,"$ARGV[1]" or die;
open TIME,">$ARGV[2]" or die;
print TIME "Element\tDistance\tInsertion_time\n";

while(<NEW>){
	chomp;
	my @array = split /\s+/;
	if($array[0] ne $array[1]){
		$hash_cal{$array[0]}++;
		if($hash_cal{$array[0]} == 1){
			open NEW1,">./$array[0].fasta" or die;
			########  reverse strand
			if($array[-4] > $array[-3]){
				my $seq = reverse_complement($hash_sequence{$array[1]});
				print NEW1 ">$array[0]\n$hash_sequence{$array[0]}\n>$array[1]\n$seq\n";
			}
			else{
				print NEW1 ">$array[0]\n$hash_sequence{$array[0]}\n>$array[1]\n$hash_sequence{$array[1]}\n";
			}
			close NEW1;

			system("muscle -in ./$array[0].fasta -out ./$array[0].muscle -clwstrict");
			system("distmat -sequence ./$array[0].muscle -nucmethod 2 -outfile ./$array[0].dist");


			my $dist;			
			open NEW1,"./$array[0].dist" or die;
			while(<NEW1>){
				chomp;
				next if (/Distance/ or /--/ or /Using/ or /^Gap/ or /^\s+$/ or /^$/);
				my @aa = split /\s+/;
				if($aa[-1] == 1){
					$dist = $aa[-3]/100;
					my $time = $dist/(2*1.3*10**-8);
					print TIME "$array[0]\t$dist\t$time\n";
				}
			}
			close NEW1;

			if($dist < 0.5){
				system("rm ./$array[0].muscle");
				system("rm ./$array[0].dist");
				system("rm ./$array[0].fasta");
			}
					
		}
	}
}
close NEW;
close TIME;


sub usage{
	my $die =<<DIE;
	usage : perl *.pl DNA_Helitron.fasta DNA_Helitron.fasta.blastn.DNA_Helitron.fasta.m8 output_time_table
DIE
}

sub reverse_complement{
	my ($seq) = @_;
	my $a = reverse $seq;
	$a =~ tr/AGCTagct/TCGAtcga/;
	return $a;
}
