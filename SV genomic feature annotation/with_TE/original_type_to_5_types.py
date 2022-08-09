import os
import sys


bed_file= sys.argv[1] #orignal TE type
out_bed = sys.argv[2]  # output 5 TE types

data = [line.rstrip().split("\t") for line in open(bed_file).readlines()]
out = open(out_bed, "w")
for d in data:
    if d[3] in ["DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT"]:
        out.write("\t".join(d[0:3]) + "\tDNA\n")
    elif d[3] == "DNA/Helitron":
        out.write("\t".join(d[0:3]) + "\tHelitron\n")
    elif d[3] in ["MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT"]:
        out.write("\t".join(d[0:3]) + "\tMITE\n")
    elif d[3] == "LTR/Copia":
        out.write("\t".join(d[0:3]) + "\tCopia\n")
    elif d[3] == "LTR/Gypsy":
        out.write("\t".join(d[0:3]) + "\tGypsy\n")
    else:
        out.write("\t".join(d[0:3]) + "\tOthers\n")

