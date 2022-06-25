import os
import sys

gene_express = "seedling.normalized.expression_sorted"
gene_coords = {}

data = [line.rstrip().split("\t") for line in open("mo17_genes_4cols.txt", "r").readlines()]

for d in data:
    gene_coords[d[3]] = d[0]+ "\t" +d[1] +"\t" +d[2]

out = open("phenotypes.bed", "w")

data = [line.split("\t") for line in open(gene_express, "r").readlines()]

out.write("#Chr\tstart\tend\tID\t" + "\t".join(data[0][1:]))
for d in data[1:]:
    out.write(gene_coords[d[0]] + "\t" + "\t".join(d))
out.close()
