#!/bin/bash
###Make sure genomic features of Mo17 genome were not overlapped.
bedtools intersect -a ../maize_up5K_uniq  -b  ../maize_up5K_uniq -wao  |awk '$2 != $5' -
bedtools intersect -a ../maize_exon_uniq  -b  ../maize_exon_uniq -wao  |awk '$2 != $5' -
bedtools intersect -a ../maize_intron_uniq -b  ../maize_intron_uniq -wao  |awk '$2 != $5' -
bedtools intersect -a ../maize_down5K_uniq -b  ../maize_down5K_uniq -wao  |awk '$2 != $5' -
bedtools intersect -a ../maize_intergenic_uniq -b  ../maize_intergenic_uniq -wao  |awk '$2 != $5' -
bedtools intersect -a ../maize_TE_uniq         -b  ../maize_TE_uniq -wao  |awk '$2 != $5' -


###Convert the postions of DEL and INS into bed files
perl convert_VCF_to_bed_file_for_INS2.pl    INS.SViper.out.1000.vcf.uniqe.vcf INS_position_info
perl convert_VCF_to_bed_file_v2_for_DEL2.pl DEL.SViper.out.1000.vcf.uniqe.vcf DEL_position_info


###Check the overlap between Mo17 genomic features and DEL & INS:
 intersectBed -a INS_position_info -b ../maize_up5K_uniq             -wo >INS_gene_up5K
 intersectBed -a INS_position_info -b ../maize_exon_uniq             -wo >INS_gene_exon
 intersectBed -a INS_position_info -b ../maize_intron_uniq           -wo >INS_gene_intron
 intersectBed -a INS_position_info -b ../maize_down5K_uniq           -wo >INS_gene_down5K
 intersectBed -a INS_position_info -b ../maize_intergenic_uniq       -wo >INS_intergenic 
 intersectBed -a INS_position_info -b ../maize_TE_uniq               -wo >INS_TE

 intersectBed -a DEL_position_info -b ../maize_up5K_uniq             -wo >DEL_gene_up5K
 intersectBed -a DEL_position_info -b ../maize_exon_uniq             -wo >DEL_gene_exon
 intersectBed -a DEL_position_info -b ../maize_intron_uniq           -wo >DEL_gene_intron
 intersectBed -a DEL_position_info -b ../maize_down5K_uniq           -wo >DEL_gene_down5K
 intersectBed -a DEL_position_info -b ../maize_intergenic_uniq       -wo >DEL_intergenic 
 intersectBed -a DEL_position_info -b ../maize_TE_uniq               -wo >DEL_TE

###The location of SVs was determined as the most overlapped genomic feature performed by bedtools.
perl get_SV_distribution2.pl INS_TE INS_gene_exon INS_gene_intron INS_gene_up5K  INS_gene_down5K INS_intergenic INS
perl get_SV_distribution2.pl DEL_TE DEL_gene_exon DEL_gene_intron DEL_gene_up5K  DEL_gene_down5K DEL_intergenic DEL


