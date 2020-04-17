import os

folder = "/mnt/c/Users/Beibei Xin/OneDrive - University of Southern California/Documents/Projects/nanoporesv/data/bxin/"

query = [line.rstrip().split("\t")[0] for line in open(folder + "phylogenetic/SViper_106.mdist.id").readlines()]

library = [line.split("\t") for line in open(folder + "106_materials_full_info.txt").readlines()]

out = open(folder+"phylogenetic/sviper_name_mapping.id","w")
name_mapping={}
type_mapping={}
for each_line in library:
	name_mapping[each_line[0]] = each_line[1]
	type_mapping[each_line[0]] = each_line[2]

for q in query:
	if q in name_mapping.keys():
		out.write(q + "\t" + name_mapping[q] + "\n")
out.close()