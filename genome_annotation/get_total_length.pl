use strict;

my($file)=@ARGV;

open F,$file or die;
my $total;
while(<F>){
	chomp;
	my @arr=split;
	my $length=$arr[2]-$arr[1]+1;
	if($length <= 0){
		die;
	}
	$total+=$length;
}
print $file ,"\t",$total,"\n";

