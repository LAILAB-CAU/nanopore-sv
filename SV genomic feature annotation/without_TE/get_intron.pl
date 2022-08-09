#!/usr/bin/perl

open I, "$ARGV[0]";

my $i=0;
while(<I>){
	chomp;
	@ss=split;
	$i++;
	$hash{$i}=$_;
}


for ($j=1;$j<=$i;$j++){
	unless(exists $hash{$j+1}){
		next;
	}
	@hh=split /\t/,$hash{$j};
	@mm=split /\t/,$hash{$j+1};
	if($hh[2]<$mm[1] && $hh[3] eq $mm[3]){
		$hh[2]=$hh[2]+1;
		$mm[1]=$mm[1]-1;
		print "$hh[0]\t$hh[2]\t$mm[1]\t$hh[3]\t$hh[4]\n";
	}
}

