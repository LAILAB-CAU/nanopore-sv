#!/bin/bash

## prepare all exon regions
data_folder=/public1/home/sc30797/bxin/01_nanopore/data/mo17_gene_exon_intron_up2k/
awk '$3=="exon"{print $_}' /public1/home/sc30797/david/code/sv/sample/Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 > $data_folder/mo17_exons.txt
grep UTR /public1/home/sc30797/david/code/sv/sample/Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 > $data_folder/mo17_UTRs.txt

## prepare mo17_genes_4cols.txt, mo17_genes_downstream_5kb.txt,  mo17_genes_upstream_5kb.txt, mo17_exons_4cols.txt,  mo17_UTRs_4cols_no_transcript.txt, mo17_UTRs_4cols.txt
python convert_to_4cols.py

# prepare mo17_introns_4cols.txt
python get_intron_regions_from_gff3.py  


#### extract UTRs from exons regions
bedtools intersect -wa -wb -a ../data/mo17_gene_exon_intron_up2k/mo17_exons_4cols.txt -b ../data/mo17_gene_exon_intron_up2k/mo17_UTRs_4cols.txt | awk '$4==$8{print $_}' > ../data/matrixeQTL/temp  # make sure transcript are unique

# prepare exons_exclude_UTRs.txt
python exclude_UTRs.py  

bedtools intersect -wa -v -a ../data/mo17_gene_exon_intron_up2k/mo17_exons_4cols.txt -b ../data/mo17_gene_exon_intron_up2k/mo17_UTRs_4cols.txt | awk -F"_" '{print $1"\texon"}' - | cat - $data_folder/exons_exclude_UTRs.txt > $data_folder/mo17_exons_4cols_no_transcript.txt


## add 4th column to exon, intron, up5kb, down5kb, UTR table
awk '{print $_}' $data_folder/mo17_exons_4cols_no_transcript.txt | sortBed -i - | uniq - > ../data/matrixeQTL/mo17_exon_uniq.txt
awk '{print $_"\tUTR"}' $data_folder/mo17_UTRs_4cols_no_transcript.txt | sortBed -i - | uniq - > ../data/matrixeQTL/mo17_UTRs_uniq.txt
awk '{print $_"\tintron"}' $data_folder/mo17_introns_4cols.txt | sortBed -i - | uniq - > ../data/matrixeQTL/mo17_intron_uniq.txt
awk '{print $_"\tup2k"}' $data_folder/mo17_genes_upstream_5kb.txt | sortBed -i - | uniq - > ../data/matrixeQTL/mo17_up5k_uniq.txt
awk '{print $_"\tdown2k"}' $data_folder/mo17_genes_downstream_5kb.txt | sortBed -i - | uniq - > ../data/matrixeQTL/mo17_down5k_uniq.txt

