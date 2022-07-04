import os
import subprocess
import os.path
import numpy as np
from operator import itemgetter
from functools import reduce
from operator import add


data = [line.rstrip().split("\t") for line in open("../data/mo17_gene_exon_intron_up2k/mo17_exons_4cols.txt").readlines()]
exons = {}
chrs = {}

out = open("../data/mo17_gene_exon_intron_up2k/mo17_introns_4cols.txt", "w")
for d in data:
    if d[3] in exons:  # if the key already exist
            exons[d[3]] += [[int(d[1]), int(d[2])]]  # form a 2D list
    else:
        exons[d[3]] = [[int(d[1]), int(d[2])]]
        chrs[d[3]] = d[0]

for gene in exons.keys():
    exons[gene] = reduce(add, sorted(exons[gene], key=itemgetter(0)))   # sort the 2D list through id 0, and reduce them to a 1D list
    if len(exons[gene]) == 2:
        continue
    else:
        pairs = int(len(exons[gene])/2)
        if len(exons[gene])%2 !=0:
            print (gene)
        for i in np.arange(pairs-1):
            out.write(chrs[gene] + "\t" + str(exons[gene][2*i+1] + 1) + "\t" + str(exons[gene][2*i+2] -1)+ "\t" + gene.split("_")[0] + "\n")
out.close()

