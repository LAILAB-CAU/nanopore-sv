import os
import sys

gvcf = sys.argv[1]

# 1       31227906        .       A       <NON_REF>       .       .       ASM_Chr=1;ASM_End=29891785;ASM_Start=29891785;ASM_Strand=+;END=31227906 GT:AD:DP:PL     0:30,0:30:0,90,90

data = [line.split("\t") for line in open(gvcf).readlines()]

del_out = open(gvcf[:-4] + "del", "w")
ins_out = open(gvcf[:-4] + "ins", "w")

for d in data:
    if len(d[3]) >= 50:
        del_out.write("chr" + d[0] + "\t" + d[1] + "\t" + str(int(d[1])+ len(d[3])) + "\t" + str(len(d[3])) + "\t"  + d[3] + "\n")
    else:
        ins_seq = d[4].split(",")[0]
        ins_out.write("chr" + d[0] + "\t" + d[1] + "\t" + str(int(d[1])+1) + "\t" + str(len(ins_seq)) + "\t"  + ins_seq + "\n")
del_out.close()
ins_out.close()