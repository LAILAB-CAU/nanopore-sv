die if ($#ARGV!=1);
open (F, $ARGV[0]);
open (O1, ">$ARGV[1].DEL.hom.vcf");
open (O2, ">$ARGV[1].INS.hom.vcf");

while (<F>) {
	chomp;
	# $fh="O$i";
	# for ($i=1;$i<=12;$i++) {
		# print $fh "$_\n" if(/^#/);
	# }
	next if(/^#/);
	if (/SVTYPE=(\S+?);/) {
		$type=$1;
	}
	if ($type eq "TRA"){
		$bnd++;
		next;
	}
	if (/SVLEN=(\S+)/) {
		$len=abs($1);
#		next if($len<50);
	}
	$total++;
	@F=split /\t/, $_;
	if ($F[6] eq "UNRESOLVED"){
		$failed++;
		next;
	}
	if (/IMPRECISE;/) {
		$imprecise++;
#		next;
	 }
	if($type eq "DEL/INV"){
		$del_inv++;
		next;
	}
	if($type eq "DUP"){
		$dup++;
		next;
	}
	if($type eq "DUP/INS"){
		$dup_ins++;
		next;
	}
	if($type eq "INV"){
		$inv++;
		next;
	}
	if($type eq "INVDUP"){
		$invdup++;
		next;
	}
	if($type eq "INV/INVDUP"){
		$dup_invdup++;
		next;
	}
	$F[8]="GT:DR:DV";
	($genotype,$ref,$alt)=split /:/, $F[9];
	if($type eq "DEL"){
		if ($genotype eq "0/1") {
			$delhet++;
		}
		elsif ($genotype eq "1/1") {
			if ($ref==0 && $alt>=5) {
				$delhomq++;
#				print O1 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			elsif ($ref==0 && $alt<5) {
				$delhomuq1++;
#				print O1 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			else{
				$delhomuq2++;
#				print O1 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
		}
		elsif ($genotype eq "0/0") {
			$delhom0++;
		}
		else{
		}
		next;
	}
	if($type eq "INS"){
		$ins++;
		if ($genotype eq "0/1") {
			$inshet++;
		}
		elsif ($genotype eq "1/1") {
			$inshom++;
			if ($ref==0 && $alt>=5) {
				$inshomq++;
#				print O2 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			elsif ($ref==0 && $alt<5) {
				$inshomuq1++;
#				print O2 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			else{
				$inshomuq2++;
#				print O2 "$F[0]\t$F[1]\t$type\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
		}
		elsif ($genotype eq "0/0") {
			$inshom0++;
		}
		else{
		}
		next;
	}
}
print "$ARGV[0]\t$total\t$failed\t$imprecise\t$bnd\t$del_inv\t$dup\t$dup_ins\t$inv\t$invdup\t$dup_invdup\t$delhet\t$delhom0\t$delhomq\t$delhomuq1\t$delhomuq2\t$inshet\t$inshom0\t$inshomq\t$inshomuq1\t$inshomuq2\n";
close F;
close O1;
close O2;






