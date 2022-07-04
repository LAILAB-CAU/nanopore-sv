import os
import subprocess
import os.path

### seperate UTR and exon regions, for those exons that overlap with UTRs
def unoverlapped_region(r1, e1, r2, e2):
    if r2-r1 >=2:
        return [r1, r2]
    elif e1-e2 >=2:
        return [e2, e1]
    else:
        return -1

data = [line.rstrip().split("\t") for line in open("../data/matrixeQTL/temp").readlines()]
out = open("../data/mo17_gene_exon_intron_up2k/exons_exclude_UTRs.txt", "w")
for d in data:
    gene_name = d[3].split("_")[0]
    if unoverlapped_region(int(d[1]), int(d[2]), int(d[5]), int(d[6])) != -1:
        out.write(d[0] + "\t" + str(unoverlapped_region(int(d[1]), int(d[2]), int(d[5]), int(d[6]))[0]) + "\t"  \
        + str(unoverlapped_region(int(d[1]), int(d[2]), int(d[5]), int(d[6]))[1]) + "\t" + gene_name + "\texon" + "\n")
out.close()
