use strict;

my($file,$file_out)=@ARGV;

open F,$file or die;
open O,">",$file_out or die;

my %info;
my %length;
my $id=0;
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
	my $start=$arr[1];

	my $length=-1;
	my @temp=split /;/,$arr[7];
	if($temp[2]=~ s/SVLEN=//g){
		if($temp[2] <0){
			$temp[2]=-$temp[2];
		}
		$length=$temp[2];
	}
	else{
		print "$arr[0]\t$arr[1]\n";
		next;
	}

	$length=int($length);
	my $end=$arr[1]+$length; #for half open interval
	my $chr=$arr[0];
	if(exists $info{$id}){
		if($length >  $length{$id}){
			$length{$id}=$length;	
			$info{$id} = "$chr\t$start\t$end\t$length{$id}\n";
		}
	}
	else{
		$length{$id}=$length;
		$info{$id} = "$chr\t$start\t$end\t$length{$id}\n";
	}
}

foreach my $key(sort keys %info){
 
	print O $info{$key};
}

