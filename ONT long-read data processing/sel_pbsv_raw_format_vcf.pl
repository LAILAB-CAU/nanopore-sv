die if ($#ARGV!=1);
open (F, $ARGV[0]);
open (O1a, ">$ARGV[1].DEL.het.vcf");
open (O1b, ">$ARGV[1].DEL.homq.vcf");
open (O1c, ">$ARGV[1].DEL.homuq1.vcf");
open (O1d, ">$ARGV[1].DEL.homuq2.vcf");
open (O1e, ">$ARGV[1].DEL.other.vcf");
open (O2a, ">$ARGV[1].INS.het.vcf");
open (O2b, ">$ARGV[1].INS.homq.vcf");
open (O2c, ">$ARGV[1].INS.homuq1.vcf");
open (O2d, ">$ARGV[1].INS.homuq2.vcf");
open (O2e, ">$ARGV[1].INS.other.vcf");
open (O3, ">$ARGV[1].INV.vcf");
open (O4, ">$ARGV[1].CNV.vcf");
open (O5, ">$ARGV[1].DUP.vcf");
open (O6, ">$ARGV[1].BND.vcf");
while (<F>) {
	chomp;
	$fh="O$i";
	for ($i=1;$i<=12;$i++) {
		print $fh "$_\n" if(/^#/);
	}
	next if(/^#/);
	if (/SVLEN=(\S+)/) {
		$len=abs($1);
		next if($len<50);
	}
	$total_number++;
	@F=split /\t/, $_;
	if ($F[6] eq "Decoy" || $F[6] eq "NotFullySpanned"){
		$failed_number++;
		next;
	}
	if ($F[2]=~/BND/){
		print O6 "$_\n";
		$bnd_num++;
		next;
	}
	if($F[2]=~/CNV/){
		print O4 "$_\n";
		$cnv_num++;
		next;
	}
	if($F[2]=~/DUP/){
		print O5 "$_\n";
		$dup_num++;
		next;
	}
	if($F[2]=~/INV/){
		print O3 "$_\n";
		$inv_num++;
		next;
	}

	$F[8]="GT:DR:DV";
	@tmp=split /:/, $F[9];
	($ref,$alt)=split /,/,$tmp[1];
	$F[9]="$tmp[0]:$ref:$alt";
	
	if($F[2]=~/DEL/){
		if ($tmp[0] eq "0/1") {
			$delhet_num++;
			print O1a "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
		}
		elsif ($tmp[0] eq "1/1") {
			if ($ref==0 && $alt>=5) {
				$delhomq_num++;
				print O1b "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			elsif ($ref==0 && $alt<5) {
				$delhomuq1_num++;
				print O1c "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			else{
				$delhomuq2_num++;
				print O1d "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
		}
		else{
			$delother_num++;
			print O1e "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
		}
		next;
	}
	if($F[2]=~/INS/){
		$ins_num++;
		if ($tmp[0] eq "0/1") {
			$inshet_num++;
			print O2a "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
		}
		elsif ($tmp[0] eq "1/1") {
			$inshom_num++;
			if ($ref==0 && $alt>=5) {
				$inshomq_num++;
				print O2b "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			elsif ($ref==0 && $alt<5) {
				$inshomuq1_num++;
				print O2c "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
			else{
				$inshomuq2_num++;
				print O2d "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
			}
		}
		else{
			$insother_num++;
			print O2e "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\n";
		}
		next;
	}
}
print "$ARGV[0]\t$total_number\t$failed_number\t$bnd_num\t$cnv_num\t$dup_num\t$inv_num\t$delother_num\t$delhet_num\t$delhomq_num\t$delhomuq1_num\t$delhomuq2_num\t$insother_num\t$inshet_num\t$inshomq_num\t$inshomuq1_num\t$inshomuq2_num\n";
close F;
close O1;
close O2;
close O3;
close O4;
close O5;
close O6;
close O7;
close O8;
close O9;
close O10;
close O11;
close O12;





