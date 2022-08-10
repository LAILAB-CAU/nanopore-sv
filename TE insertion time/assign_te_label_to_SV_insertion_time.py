import os
import sys
import numpy as np
import subprocess


intersect_output = sys.argv[1]
final_output = sys.argv[2]

out = open(final_output, "w")

data = [line.rstrip().split("\t") for line in open(intersect_output, "r").readlines()]

for d in data:
    if len(d) > 3:  
        out.write(d[0] + "\t" + d[1] + "\t" + d[2] + "\t")
        tes = [each.split(":")[0] for each in d[3][:-1].split("|")]
        lens = [int(each.split(":")[1]) for each in d[3][:-1].split("|")]
        insertion_time = [each.split(":")[2] for each in d[3][:-1].split("|")]
        if len(tes) == 1:   
            out.write(tes[0] + "\t" + insertion_time[0] + "\n")
        else:
            index = np.argmax(lens)
            out.write(tes[index] + "\t" + insertion_time[index] + "\n")
out.close()