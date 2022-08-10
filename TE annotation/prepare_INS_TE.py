import os
import sys
import numpy as np
import subprocess


INS_id = sys.argv[1]
# ./section4/EDTA/INS_id.bed
# 22994	chr10_1269462_1269711
# 38102	chr2_46752600_46752932
# 130031	chr8_22821289_22821360
# 65713	chr3_212685152_212685203
INS_fa_out = sys.argv[2]

# ./section4/EDTA/INS.fa.out
#    SW   perc perc perc  query     position in query    matching                   repeat                 position in repeat
# score   div. del. ins.  sequence  begin end   (left)   repeat                     class/family       begin   end    (left)      ID

#   880    4.2  0.0  5.8  1            13   140    (0) C opie_AC210610_10772        LTR/Copia           (3458)   4024    3904      1  
#   567   25.1  5.3  2.0  10           11   254    (7) C TE_00002629                MITE/DTT               (0)    252       1      2 
final_output = sys.argv[3]
# chr1    1322    1469    DNA/DTA
# chr1    1347    1527    MITE/DTT
# chr1    1470    1552    MITE/DTH

def te_type(te_string):
    te_result = ""
    if te_string in ["DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT"]:
        te_result="DNA"
    elif te_string == "DNA/Helitron":
        te_result="Helitron"
    elif te_string in ["MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT"]:
        te_result="MITE"
    elif te_string == "LTR/Copia":
        te_result="Copia"
    elif te_string == "LTR/Gypsy":
        te_result="Gypsy"
    else:
        te_result="Others"
    return te_result


out = open(final_output, "w")

data = [line.rstrip().split() for line in open(INS_id, "r").readlines()]

ins_id = {}
for d in data:
    items = d[1].split("_")
    ins_id[d[0]] = "\t".join(items)

data = [line.rstrip().split() for line in open(INS_fa_out, "r").readlines()[3:]]

for d in data:
    ins_chr, ins_start, ins_end = ins_id[d[4]].split("\t")
    out.write(ins_id[d[4]]+ "\t")  
    out.write(ins_chr + "\t" + str(int(ins_start) + int(d[5])) + "\t" + str(int(ins_start) + int(d[6])) + "\t")
    out.write(te_type(d[10]) + "\t"  + str(int(d[6])-int(d[5])+1) + "\n")
out.close()
