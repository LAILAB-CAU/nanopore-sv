import sys
import os
import numpy as np


# prepare fastQTL batch files
input_file = sys.argv[1]
k = int(sys.argv[2]) # separate the fastQTL commands to k different .sh files

out_file_prefix = ".".join(input_file.split(".")[:-1])

data = open(input_file).readlines()

if len(data)%k == 0:
    chunk_size = int(len(data)/k)
else:
    chunk_size= int(len(data)/k) +1

for i in np.arange(k):
    data_to_file = data[(i*chunk_size): ((i+1)*chunk_size)]
    out_file = out_file_prefix + ".chunk" + str(i) + ".sh"
    out = open(out_file, "w")
    out.write("#!/bin/bash\n\n")
    out.write("#SBATCH --ntasks=1 \n#SBATCH --mem=60GB\n#SBATCH --partition=amd_256\n#SBATCH --time=120:00:00\n#SBATCH --error=/public1/home/sc30797/bxin/04_Effects_of_genomic_variations/system_io/%x-%j.e\n#SBATCH --output=/public1/home/sc30797/bxin/04_Effects_of_genomic_variations/system_io/%x-%j.o\n")
    out.write("conda activate\n")
    out.write("export LD_LIBRARY_PATH=/public1/home/sc30797/bxin/softwares/caviar-master/CAVIAR-C++/:$LD_LIBRARY_PATH\n\n")
    out.write("".join(data_to_file))
    out.write("wait\n")
    out.close()
    cmd = "sbatch " + out_file
    os.system(cmd)

