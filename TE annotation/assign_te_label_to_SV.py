import os
import sys
import numpy as np
import subprocess


intersect_output = sys.argv[1]
final_output = sys.argv[2]

out = open(final_output, "w")

data = [line.rstrip().split("\t") for line in open(intersect_output, "r").readlines()]

for d in data:
    out.write(d[0] + "\t" + d[1] + "\t" + d[2] + "\t")
    if len(d) == 3:  
        out.write("non-TE\n")
    else:
        tes = [each.split(":")[0] for each in d[3][:-1].split("|")]
        lens = [int(each.split(":")[1]) for each in d[3][:-1].split("|")]
        if len(tes) == 1 and tes[0] == "Others":  
            out.write("Others\n")
        elif len(tes) == 1:   
            out.write(tes[0] + "\n")
        elif "Others" in tes:  
            index = tes.index("Others")
            tes.pop(index)
            lens.pop(index)
            index = np.argmax(lens)
            out.write(tes[index] + "\n")
        else:
            index = np.argmax(lens)
            out.write(tes[index] + "\n")
out.close()