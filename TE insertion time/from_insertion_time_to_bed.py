import sys
import os
import subprocess

input_file = sys.argv[1]  # insertion time for LTRs 
output_file = sys.argv[2]

def search_te_info(in_file, id):
    # in_file is the bedgraph file from eQTL analyses, sv_info is "chr1\t45667"
    cmd = subprocess.Popen(['grep', '-P', id, in_file], 
           stdout=subprocess.PIPE, 
           stderr=subprocess.STDOUT)
    stdout,stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    if 'grep' in output:
        print(" ".join(['grep', '-P', id, in_file]))
        print("search process error!")
        return "-1"
    else:
        return output.split("\t")


#check the p-value for each gene
data =[line.rstrip().split("\t") for line in open(input_file).readlines()[1:]]


out = open(output_file, "w")
query_file = "Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.intact.gff3"
for d in data:
    id = "ID=" + d[0].split(";")[0]
    te = d[0].split("Classification=")[1].split(";")[0]
    query_result = search_te_info(query_file, id)
    out.write(query_result[0] + "\t" + query_result[3] + "\t" + query_result[4] + "\t" + te + "\t" + d[2] + "\n")
out.close()