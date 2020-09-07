#!/usr/bin/env perl
use strict;
use warnings;

open LIST, "./bamlist" or die;


my $cal;
my $file_cal;
my $vcfref = "/public1/home/sc30797/ont_data/ngs_bams/bcf_snp/bcftools.vcf.gz";
while(my $file = <LIST>){
	chomp $file;
    my @a = split /\//,$file;
	my @b = split /\./,$a[-1];	
	my $sample = $b[0];
	$cal++;
	my $a = $cal % 16;
	if($a == 1){
		$file_cal++;
		# open NEW,">./BaseRecalibrator_$file_cal.sh" or die;
		# print NEW "#!/bin/bash\n";
		# print NEW "gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk\n";
		# print NEW "refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa\n";
		# print NEW "path=~/ont_data/ngs_bams/bqsr\n";
		# print NEW "\$gatk BaseRecalibrator \\\n -R \$refgenome \\\n -I $file \\\n --known-sites $vcfref \\\n -O \$path/$sample\.sorted.markdup.recal_data.table && echo \"** $sample\.sorted.markdup.recal_data.table done **\" &\n";

        open NEW2, ">./ApplyBQSR_$file_cal.sh" or die;
        print NEW2 "#!/bin/bash\n";
        print NEW2 "gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk\n";
		print NEW2 "refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa\n";
        print NEW2 "path=~/ont_data/ngs_bams/bqsr\n";
        print NEW2 "\$gatk ApplyBQSR \\\n --bqsr-recal-file \$path/$sample\.sorted.markdup.recal_data.table \\\n -R \$refgenome \\\n -I $file \\\n -O \$path/$sample\_marked.sorted.BQSR.bam && echo \"** ApplyBQSR done **\" &\n";
    }
	elsif($a == 0){
		# print NEW "\$gatk BaseRecalibrator \\\n -R \$refgenome \\\n -I $file \\\n --known-sites $vcfref \\\n -O \$path/$sample\.sorted.markdup.recal_data.table && echo \"** $sample\.sorted.markdup.recal_data.table done **\" &\n";
		# print NEW "wait\n";	

        print NEW2 "\$gatk ApplyBQSR \\\n --bqsr-recal-file \$path/$sample\.sorted.markdup.recal_data.table \\\n -R \$refgenome \\\n -I $file  \\\n -O \$path/$sample\_marked.sorted.BQSR.bam && echo \"** ApplyBQSR done **\" &\n";
        print NEW2 "wait\n";
        print NEW2 "echo \"** ApplyBQSR $file_cal done **\"";
    }
	else{
		# print NEW "\$gatk BaseRecalibrator \\\n -R \$refgenome \\\n -I $file \\\n --known-sites $vcfref \\\n -O \$path/$sample\.sorted.markdup.recal_data.table && echo \"** $sample\.sorted.markdup.recal_data.table done **\" &\n";
        print NEW2 "\$gatk ApplyBQSR \\\n --bqsr-recal-file \$path/$sample\.sorted.markdup.recal_data.table \\\n -R \$refgenome \\\n -I $file  \\\n -O \$path/$sample\_marked.sorted.BQSR.bam && echo \"** ApplyBQSR done **\" &\n";
        }
	
}
