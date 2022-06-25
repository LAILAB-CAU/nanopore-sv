import os
import sys

# egene_list = sys.argv[1]

# genes = [line.rstrip() for line in open(egene_list).readlines()]

# for g in genes:
# 	cmd = "rm /data05/bxin/contribution_of_variations/combined_eGenes/" + g + ".*"
# 	os.system(cmd)


data = [line.split("\t") for line in open("/phenotypes.bed").readlines()]
chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "#Chr"]

out = open("phenotypes.bed.new", "w")
for d in data:
	if d[0] in chrs:
		out.write("\t".join(d))
out.close()
