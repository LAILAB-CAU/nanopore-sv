import os
import sys

input_file = "/public1/home/sc30797/bxin/04_Effects_of_genomic_variations/SV_SNP_INDEL/combined.all.chunks.CAVIAR.results.txt"
out = open("/public1/home/sc30797/bxin/04_Effects_of_genomic_variations/SV_SNP_INDEL/combined.all.chunks.CAVIAR.results.2cols.txt","w")

data = [line.rstrip().split("\t") for line in open(input_file).readlines()]

for d in data:
    if len(d) == 2:
        temp_gene_id = d[0]
        out.write(d[0] + "\t" + d[1] + "\n")
    else:
        out.write(temp_gene_id + "\t" + d[0] + "\n")
out.close()
