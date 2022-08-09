#!/bin/bash

gff=Zm-Mo17-REFERENCE-CAU-1.0_mgdb.gff3
edta_result=mo17.EDTA.5.TE.type.bed

##Obtain the genomic positions of exon,intron,gene,Upstream5K,Downstream5K,TE.
perl get_gene_bed.pl $gff  Mo17_CAU_gene_bed
perl get_exon_bed_v2.pl $gff  Mo17_CAU_exon_bed

cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}if($ss[4] eq "-"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}' >maize_Mo17_up5K.bed

cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}if($ss[4] eq "-"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}' >maize_Mo17_down5K.bed

cat Mo17_CAU_exon_bed  | sort -k4,4 -k1,1 -k2,2n > maize_Mo17_exon.bed.sorted
cat maize_Mo17_exon.bed.sorted | awk '{print $1":"$4"\t"$2"\t"$3}' | bedtools merge -i - > maize_Mo17_exon.bed.tmp
cat maize_Mo17_exon.bed.tmp|perl -ne 'chomp;@ss=split;@mm=split /:/,$ss[0];print "$mm[0]\t$ss[1]\t$ss[2]\t$mm[1]\n";' >maize_Mo17_exon.bed.merge
perl get_intron.pl maize_Mo17_exon.bed.merge >maize_Mo17_intron.bed

cat  $edta_result | awk '{print $1"\t"$2"\t"$3}'|grep \# -v|sort -k1,1 -k2,2n | bedtools merge -i - > maize_Mo17_TE_merge


##Split the positions of TE，exon ,intron ,Upstream5K,Downstream5K into chromosome-level

for((i=1;i<=10;i++));  
do  
mkdir chr$i
perl get_chr.pl maize_Mo17_down5K.bed $i  > chr$i/maize_Mo17_down5K.bed_chr$i
perl get_chr.pl maize_Mo17_intron.bed $i  > chr$i/maize_Mo17_intron.bed_chr$i
perl get_chr.pl maize_Mo17_up5K.bed   $i  > chr$i/maize_Mo17_up5K.bed_chr$i
perl get_chr.pl Mo17_CAU_exon_bed     $i  > chr$i/Mo17_CAU_exon_bed_chr$i
perl get_chr.pl maize_Mo17_TE_merge   $i  > chr$i/Mo17_CAU_TE_bed_chr$i
perl get_chr.pl chrome_length         $i  > chr$i/chrome_length_chr$i
done


##Assign a label to each nucleotide according to the priority of TE, exon, intron, up5K, down5K, intergenic.

for((i=1;i<=10;i++));  
do  
perl get_interval.pl chr$i/Mo17_CAU_TE_bed_chr$i chr$i/Mo17_CAU_exon_bed_chr$i chr$i/maize_Mo17_intron.bed_chr$i chr$i/maize_Mo17_up5K.bed_chr$i chr$i/maize_Mo17_down5K.bed_chr$i  chr$i/chrome_length_chr$i  >  chr$i/interval & 
done

wait

for((i=1;i<=10;i++));  
do  

###Merge adjacent nucleotides with the same genomic feature.
cat chr$i/interval|grep TE        | perl merge_file.pl -      >chr$i/maize_TE_uniq & 
cat chr$i/interval|grep exon      | perl merge_file.pl -      >chr$i/maize_exon_uniq &
cat chr$i/interval|grep intron    | perl merge_file.pl -      >chr$i/maize_intron_uniq &
cat chr$i/interval|grep up2K      | perl merge_file.pl -      >chr$i/maize_up5K_uniq &
cat chr$i/interval|grep down2K    | perl merge_file.pl -      >chr$i/maize_down5K_uniq &

cat chr$i/interval|grep intergenic| perl merge_file.pl -      >chr$i/maize_intergenic_uniq &

wait
done

###Merge results from different chromosomes.
cat chr*/maize_TE_uniq       > maize_TE_uniq
cat chr*/maize_exon_uniq       > maize_exon_uniq
cat chr*/maize_intron_uniq     > maize_intron_uniq
cat chr*/maize_up5K_uniq       > maize_up5K_uniq
cat chr*/maize_down5K_uniq     > maize_down5K_uniq
cat chr*/maize_intergenic_uniq > maize_intergenic_uniq


