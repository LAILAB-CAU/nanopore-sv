import os
import subprocess
import os.path

data = [line.rstrip().split("\t") for line in open("../data/matrixeQTL/mo17_genes.txt").readlines()]

out = open("../data/matrixeQTL/mo17_genes_upstream_5kb.txt", "w")
for d in data:
    gene_name = d[8].split(";")[0].split("=")[1]
    if d[6] == "+":
        if (int(d[3])-2000) >0:
            out.write(d[0] + "\t" + str(int(d[3])-5000) + "\t" + d[3] + "\t" + gene_name + "\n")
    else:
        out.write(d[0] + "\t" + d[4] + "\t" + str(int(d[4])+5000) + "\t" + gene_name + "\n")
out.close()

out = open("../data/mo17_gene_exon_intron_up2k/mo17_genes_downstream_5kb.txt", "w")
for d in data:
    gene_name = d[8].split(";")[0].split("=")[1]
    if d[6] == "-":
        if (int(d[3])-2000) >0:
            out.write(d[0] + "\t" + str(int(d[3])-5000) + "\t" + d[3] + "\t" + gene_name + "\n")
    else:
        out.write(d[0] + "\t" + d[4] + "\t" + str(int(d[4])+5000) + "\t" + gene_name + "\n")
out.close()

out = open("../data/matrixeQTL/mo17_genes_4cols.txt", "w")
for d in data:
    gene_name = d[8].split(";")[0].split("=")[1]
    out.write(d[0] + "\t" + d[3] + "\t" + d[4] + "\t" + gene_name + "\n")
out.close()


data = [line.rstrip().split("\t") for line in open("../data/mo17_gene_exon_intron_up2k/mo17_exons.txt").readlines()]
out = open("../data/mo17_gene_exon_intron_up2k/mo17_exons_4cols.txt", "w")
for d in data:
    gene_name = d[8].split(";")[1].split("=")[1]
    out.write(d[0] + "\t" + d[3] + "\t" + d[4] + "\t" + gene_name + "\n")
out.close()


data = [line.rstrip().split("\t") for line in open("../data/mo17_gene_exon_intron_up2k/mo17_UTRs.txt").readlines()]
out = open("../data/mo17_gene_exon_intron_up2k/mo17_UTRs_4cols_no_transcript.txt", "w")
for d in data:
    gene_name = d[8].split("=")[1]
    out.write(d[0] + "\t" + d[3] + "\t" + d[4] + "\t" + gene_name.split("_")[0] + "\n")
out.close()

data = [line.rstrip().split("\t") for line in open("../data/mo17_gene_exon_intron_up2k/mo17_UTRs.txt").readlines()]
out = open("../data/mo17_gene_exon_intron_up2k/mo17_UTRs_4cols.txt", "w")
for d in data:
    gene_name = d[8].split("=")[1]
    out.write(d[0] + "\t" + d[3] + "\t" + d[4] + "\t" + gene_name + "\n")
out.close()
