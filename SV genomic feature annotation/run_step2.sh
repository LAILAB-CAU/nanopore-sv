perl convert_VCF_to_bed_file_for_INS2.pl    INS.vcf INS_position_info
perl convert_VCF_to_bed_file_v2_for_DEL2.pl DEL.vcf DEL_position_info

intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_up2K_uniq             -wo >INS_gene_up2K
intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_exon_uniq             -wo >INS_gene_exon
intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_intron_uniq           -wo >INS_gene_intron
intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_down2K_uniq           -wo >INS_gene_down2K
intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_intergenic_uniq       -wo >INS_intergenic 
intersectBed -a INS_position_info -b genome_different_region_for_Mo17/maize_Mo17_TE_merge          -wo >INS_TE

intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_up2K_uniq             -wo >DEL_gene_up2K
intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_exon_uniq             -wo >DEL_gene_exon
intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_intron_uniq           -wo >DEL_gene_intron
intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_down2K_uniq           -wo >DEL_gene_down2K
intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_intergenic_uniq       -wo >DEL_intergenic 
intersectBed -a DEL_position_info -b genome_different_region_for_Mo17/maize_Mo17_TE_merge         -wo >DEL_TE

perl get_SV_distribution2.pl INS_TE INS_gene_exon INS_gene_intron INS_gene_up2K  INS_gene_down2K INS_intergenic INS
perl get_SV_distribution2.pl DEL_TE DEL_gene_exon DEL_gene_intron DEL_gene_up2K  DEL_gene_down2K DEL_intergenic DEL
