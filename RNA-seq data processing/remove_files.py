import os
import sys

data = [line.split("\t") for line in open("/phenotypes.bed").readlines()]
chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "#Chr"]

out = open("phenotypes.bed.new", "w")
for d in data:
	if d[0] in chrs:
		out.write("\t".join(d))
out.close()
