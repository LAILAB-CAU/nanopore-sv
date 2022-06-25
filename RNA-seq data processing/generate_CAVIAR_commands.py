
import os
import sys
import subprocess
import numpy as np

# prepare CAVIAR batch bash files
input_file = sys.argv[1]  # full list of egenes
k = int(sys.argv[2]) # separate the GCTA egenes into k different set of genes

out_file_prefix = ".".join(input_file.split(".")[:-2])

data = open(input_file).readlines()

if len(data)%k == 0:
    chunk_size = int(len(data)/k)
else:
    chunk_size= int(len(data)/k) +1

out_file = out_file_prefix + ".sh"
all_out = open(out_file, "w")

for i in np.arange(k):
    data_to_file = data[(i*chunk_size): ((i+1)*chunk_size)]
    out_file = out_file_prefix + ".chunk" + str(i) + ".txt"
    out = open(out_file, "w")
    out.write("".join(data_to_file))
    out.close()
    all_out.write("python run_CAVIAR_for_every_eGene.py " + out_file_prefix + ".chunk" + str(i) + ".txt &\n")
all_out.close()

