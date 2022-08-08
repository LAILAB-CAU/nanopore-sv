#!/bin/bash

perl get_gene_bed.pl Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3  Mo17_CAU_gene_bed
perl get_exon_bed_v2.pl Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3  Mo17_CAU_exon_bed

cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}if($ss[4] eq "-"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}' >maize_Mo17_up5K.bed

cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}if($ss[4] eq "-"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}' >maize_Mo17_down5K.bed

cat Mo17_CAU_exon_bed  | sort -k4,4 -k1,1 -k2,2n > maize_Mo17_exon.bed.sorted
cat maize_Mo17_exon.bed.sorted | awk '{print $1":"$4"\t"$2"\t"$3}' | bedtools merge -i - > maize_Mo17_exon.bed.tmp
cat maize_Mo17_exon.bed.tmp|perl -ne 'chomp;@ss=split;@mm=split /:/,$ss[0];print "$mm[0]\t$ss[1]\t$ss[2]\t$mm[1]\n";' >maize_Mo17_exon.bed.merge
perl get_intron.pl maize_Mo17_exon.bed.merge >maize_Mo17_intron.bed

cat  ../TE annotation/mo17_EDTA_5_TE_types.txt | awk '{print $1"\t"$2"\t"$3}'|grep \# -v|sort -k1,1 -k2,2n | bedtools merge -i - > maize_Mo17_TE_merge
perl get_no_TE_region.pl maize_Mo17_TE_merge Mo17_chrom_length  >  no_TE_region_temp
cat no_TE_region_temp | perl merge_file.pl - > no_TE_region
perl get_interval.pl Mo17_CAU_exon_bed  maize_Mo17_intron.bed maize_Mo17_up5K.bed maize_Mo17_down5K.bed no_TE_region >interval

cat interval|grep exon| perl merge_file.pl - >maize_exon_uniq

cat interval|grep intron|perl merge_file.pl - >maize_intron_uniq

cat interval|grep up5K| perl merge_file.pl -  >maize_up5K_uniq

cat interval|grep down5K| perl merge_file.pl -  >maize_down5K_uniq

cat interval|grep intergenic| perl merge_file.pl - >maize_intergenic_uniq
