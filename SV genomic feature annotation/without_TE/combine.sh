#!/bin/bash
##Get the positions of exon,intron,CDS,5TUR，3UTR,gene,Upstream5K,Downstream5K
gff3=Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3

perl get_gene_bed.pl    gff3 Mo17_CAU_gene_bed
perl get_exon_bed_v2.pl gff3 Mo17_CAU_exon_bed
perl get_CDS.pl         gff3 Mo17_CAU_CDS_bed
perl get_3UTR.pl        gff3 Mo17_CAU_3UTR_bed
perl get_5UTR.pl        gff3 Mo17_CAU_5UTR_bed

cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}if($ss[4] eq "-"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}' >maize_Mo17_up5K.bed
cat Mo17_CAU_gene_bed  |perl -ne 'chomp;@ss=split;if($ss[4] eq "+"){$pos=$ss[2]+5000;$pos2=$ss[2]+1;print "$ss[0]\t$pos2\t$pos\t$ss[3]\t$ss[4]\n";}if($ss[4] eq "-"){$pos=$ss[1]-5000;$pos2=$ss[1]-1;if($pos>0){print "$ss[0]\t$pos\t$pos2\t$ss[3]\t$ss[4]\n";}else{print "$ss[0]\t0\t$pos2\t$ss[3]\t$ss[4]\n";}}' >maize_Mo17_down5K.bed

cat Mo17_CAU_exon_bed  | sort -k4,4 -k1,1 -k2,2n > maize_Mo17_exon.bed.sorted
cat maize_Mo17_exon.bed.sorted | awk '{print $1":"$4"\t"$2"\t"$3}' | /public1/home/sc30797/liuyuxin/software/bedtools2/bin/bedtools merge -i - > maize_Mo17_exon.bed.tmp
cat maize_Mo17_exon.bed.tmp|perl -ne 'chomp;@ss=split;@mm=split /:/,$ss[0];print "$mm[0]\t$ss[1]\t$ss[2]\t$mm[1]\n";' >maize_Mo17_exon.bed.merge
perl get_intron.pl maize_Mo17_exon.bed.merge >maize_Mo17_intron.bed


##Split the positions of CDS, 5TUR，3UTR,Upstream5K,Downstream5K,intron into chromosome-level
for((i=1;i<=10;i++));  
do  
perl get_chr.pl maize_Mo17_down5K.bed $i  > chr$i/maize_Mo17_down5K.bed_chr$i
perl get_chr.pl maize_Mo17_intron.bed $i  > chr$i/maize_Mo17_intron.bed_chr$i
perl get_chr.pl maize_Mo17_up5K.bed   $i  > chr$i/maize_Mo17_up5K.bed_chr$i
perl get_chr.pl Mo17_CAU_3UTR_bed     $i  > chr$i/Mo17_CAU_3UTR_bed_chr$i
perl get_chr.pl Mo17_CAU_5UTR_bed     $i  > chr$i/Mo17_CAU_5UTR_bed_chr$i
perl get_chr.pl Mo17_CAU_CDS_bed      $i  > chr$i/Mo17_CAU_CDS_bed_chr$i
perl get_chr.pl chrome_length         $i  > chr$i/chrome_length_chr$i
done

##Assign a label to every nucleotide of the genome in the priority of CDS, UTR, intron, up5K, down5K, intergenic
for((i=1;i<=10;i++));  
do  
perl get_interval.pl chr$i/Mo17_CAU_CDS_bed_chr$i chr$i/Mo17_CAU_5UTR_bed_chr$i chr$i/Mo17_CAU_3UTR_bed_chr$i chr$i/maize_Mo17_intron.bed_chr$i chr$i/maize_Mo17_up5K.bed_chr$i chr$i/maize_Mo17_down5K.bed_chr$i chr$i/chrome_length_chr$i > chr$i/interval & 
done

###Merge adjacent nucleotides with the same genomic feature
for((i=1;i<=10;i++));  
do  

cat chr$i/interval|grep 5utr      | perl merge_file.pl -      >chr$i/maize_5utr_uniq & 
cat chr$i/interval|grep 3utr      | perl merge_file.pl -      >chr$i/maize_3utr_uniq &
cat chr$i/interval|grep cds       | perl merge_file.pl -      >chr$i/maize_cds_uniq &
cat chr$i/interval|grep intron    | perl merge_file.pl -      >chr$i/maize_intron_uniq &
cat chr$i/interval|grep up5K      | perl merge_file.pl -      >chr$i/maize_up5K_uniq &
cat chr$i/interval|grep down5K    | perl merge_file.pl -      >chr$i/maize_down5K_uniq &
cat chr$i/interval|grep intergenic| perl merge_file.pl -      >chr$i/maize_intergenic_uniq &

wait
done


###Merge results from different chromosomes
cat chr*/maize_5utr_uniq       > maize_5utr_uniq
cat chr*/maize_3utr_uniq       > maize_3utr_uniq
cat chr*/maize_cds_uniq        > maize_cds_uniq
cat chr*/maize_intron_uniq     > maize_intron_uniq
cat chr*/maize_up5K_uniq       > maize_up5K_uniq
cat chr*/maize_down5K_uniq     > maize_down5K_uniq
cat chr*/maize_intergenic_uniq > maize_intergenic_uniq

