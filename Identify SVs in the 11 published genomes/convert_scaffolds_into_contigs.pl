#!/usr/bin/perl
use strict;
use warnings;

die usage() if @ARGV == 0;

my %hash_seq;
my $seq_name;
open NEW,"$ARGV[0]" or die;
while(<NEW>){
	chomp;
	if(/>/){
		$seq_name = $_;
		$seq_name =~ s/>//;
	}
	else{
		$hash_seq{$seq_name} .= $_;
	}
}
close NEW;

my $contig_count;
open NEW1,">$ARGV[1]" or die;
foreach my $key(keys %hash_seq){
	my @aa = split /N+/i,$hash_seq{$key};
	foreach my $line(@aa){
		$contig_count++;
		print NEW1 ">Contig$contig_count\n";
		print NEW1 "$line\n";
	}
}
close NEW;


sub usage{
	my $die =<<DIE;
	usage : perl *.pl scaffolds.fasta contigs.fasta
DIE
}
