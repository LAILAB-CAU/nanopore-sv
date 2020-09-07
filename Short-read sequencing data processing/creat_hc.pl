#!/usr/bin/env perl
use strict;
use warnings;

open LIST, "./bamlist" or die;


my $cal;
my $file_cal;
my @jubing = ("CHR1","CHR2","CHR3","CHR4","CHR5","CHR6","CHR7","CHR8","CHR9","CHR10");
my %handles;
while(my $file = <LIST>){
	chomp $file;
    my @a = split /\//,$file;
	my @b = split /\./,$a[-1];	
	my $sample = $b[0];
	$cal++;
	my $a = $cal % 36;
	if($a == 1){
		$file_cal++;
        for(my $i = 1;$i <=10;$i++){
            open my $jubing, ">chr$i\_$file_cal\_HaplotypeCaller.sh" or die;
            $handles{"$i-$file_cal"} =  $jubing;
            print {$handles{"$i-$file_cal"}} "#!/bin/bash\n";
            print {$handles{"$i-$file_cal"}} "gatk=/public1/home/sc30797/wangzijian/software/gatk-4.1.0.0/gatk\n";
            print {$handles{"$i-$file_cal"}} "refgenome=/public1/home/sc30797/shijunpeng/data/Mo17/Zm-Mo17-REFERENCE-CAU-1.0.fa\n";
            print {$handles{"$i-$file_cal"}} "outdir=~/ont_data/ngs_bams/gatk\n";
            print {$handles{"$i-$file_cal"}} "\$gatk  --java-options \"-Xmx4g\" HaplotypeCaller \\\n  --emit-ref-confidence GVCF \\\n -R \$refgenome \\\n -I $file \\\n -L chr$i \\\n -O \$outdir/$sample.HC.chr$i.g.vcf.gz && echo \"** $sample.HC.$i.g.vcf.gz done **\" &\n";
        }
        

    }
	elsif($a == 0){
		for(my $i = 1;$i <=10;$i++){
            print {$handles{"$i-$file_cal"}} "\$gatk  --java-options \"-Xmx4g\" HaplotypeCaller \\\n  --emit-ref-confidence GVCF \\\n -R \$refgenome \\\n -I $file \\\n -L chr$i \\\n -O \$outdir/$sample.HC.chr$i.g.vcf.gz && echo \"** $sample.HC.$i.g.vcf.gz done **\" &\n";
            print {$handles{"$i-$file_cal"}} "wait\n";
            print {$handles{"$i-$file_cal"}} "echo \"** HC chr$i $file_cal done **\"";
            close $handles{"$i-$file_cal"};
        }
    }
    elsif($a > 2 and $cal == 106){
        for(my $i = 1;$i <=10;$i++){
            print {$handles{"$i-$file_cal"}} "\$gatk  --java-options \"-Xmx4g\" HaplotypeCaller \\\n  --emit-ref-confidence GVCF \\\n -R \$refgenome \\\n -I $file \\\n -L chr$i \\\n -O \$outdir/$sample.HC.chr$i.g.vcf.gz && echo \"** $sample.HC.$i.g.vcf.gz done **\" &\n";
            print {$handles{"$i-$file_cal"}} "wait\n";
            print {$handles{"$i-$file_cal"}} "echo \"** HC chr$i $file_cal done **\"";
        }
    }
	else{
        for(my $i = 1;$i <=10;$i++){
            print {$handles{"$i-$file_cal"}} "\$gatk  --java-options \"-Xmx4g\" HaplotypeCaller \\\n  --emit-ref-confidence GVCF \\\n -R \$refgenome \\\n -I $file \\\n -L chr$i \\\n -O \$outdir/$sample.HC.chr$i.g.vcf.gz && echo \"** $sample.HC.$i.g.vcf.gz done **\" &\n";
        }
    }
	
}
close LIST;