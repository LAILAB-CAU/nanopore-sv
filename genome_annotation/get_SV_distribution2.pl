use strict;
my ($TE_file,$exon_file,$intron_file,$up_file,$down_file,$inter_file,$tag)=@ARGV;

print $tag,"\t";

open T,$TE_file or die;
open E,$exon_file or die;
open IN,$intron_file or die;
open U,$up_file or die;
open D,$down_file or die;
open I,$inter_file or die;

my %total;
my %TE;
while(<T>){
	my @arr=split;
	my $id="$arr[0]\t$arr[1]";
	$TE{$id}+=$arr[-1];
	$total{$id}++;
}
my %exon;
while(<E>){
	my @arr=split;
        my $id="$arr[0]\t$arr[1]";
        $exon{$id}+=$arr[-1];
	$total{$id}++;
}

my %intron;
while(<IN>){
	my @arr=split;
        my $id="$arr[0]\t$arr[1]";
	$intron{$id}+=$arr[-1];
        $total{$id}++;
}


my %up;
while(<U>){
	my @arr=split;
	my $id="$arr[0]\t$arr[1]";
	$up{$id}+=$arr[-1];
	$total{$id}++;
}

my %down;
while(<D>){
	my @arr=split;
        my $id="$arr[0]\t$arr[1]";
        $down{$id}+=$arr[-1];
	$total{$id}++;
}

my %inter;
while(<I>){
	my @arr=split;
        my $id="$arr[0]\t$arr[1]";
	$inter{$id}+=$arr[-1];
	$total{$id}++;
}


my $TE_num=0;
my $exon_num=0;
my $intro_num=0;
my $up_num=0;
my $down_num=0;
my $inter_num=0;
my $total_num=0;
foreach my $id( keys %total){
	my @arr=("NA",$TE{$id},$exon{$id},$intron{$id},$up{$id},$down{$id},$inter{$id});
	my $value=-10;
	my $index;
	for(my $i=1;$i<@arr;$i++){
		if($value<$arr[$i]){
			$value=$arr[$i];
			$index=$i;
		}
	}

	if($index==1){
		$TE_num++;
	}
	elsif($index==2){
		$exon_num++;
	}
	elsif($index == 3){
		$intro_num++;
	}
	elsif($index ==4){
		$up_num++;
	}
	elsif($index==5){
		$down_num++;
	}
	elsif($index==6){
		$inter_num++;
	}
	else{
		print $id,"\t",$index,"\n";
	}
	$total_num++;
}

print $TE_num,"\t",$exon_num,"\t",$intro_num,"\t",$up_num,"\t",$down_num,"\t",$inter_num,"\t",$total_num,"\n";
