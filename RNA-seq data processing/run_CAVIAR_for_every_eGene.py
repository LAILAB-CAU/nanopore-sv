import os
import sys
import subprocess
import numpy as np
from scipy.stats import norm

egene_file = sys.argv[1]
egene_list= [line.rstrip() for line in open(egene_file).readlines()]

folder = "./"
### for every gene_id, find the 1000 most significant SNPs and small InDels, and the one most signiciant SV
# egene_list= [line.rstrip() for line in open(folder + "/combined.all.chunks.benjamini.txt").readlines()]

query_file = folder + "/map_covariates.allchrs.combined.txt"

genotype_dict = {"0/0":"0", "1/1":"1","./.":"NA"}
def search_gene_name(in_file, gene_name):
    # search gene_name for each gene in
    cmd = subprocess.Popen(['grep', gene_name, in_file],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT)
    stdout,stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    if output == '':
        print(gene_name + " not found")
        return -1
    else:
        return output.rstrip().split("\n")

def seperate_two_lists(whole_list, match):
    match_list = []
    unmatch_list = []
    for word in whole_list:
        i = 0
        for match_i in match:
            if match_i in word:
                match_list.append(word)
                i = 1
        if i == 0:
            unmatch_list.append(word)
    return [item.split(" ") for item in match_list], [item.split(" ") for item in unmatch_list]

def rank_list(whole_list):
    whole_list.sort(key=lambda row: float(row[3]))

def search_sv_info(in_file, sv_id):
    cmd = subprocess.Popen(['zgrep', "-w", sv_id + "\|#CHROM" , in_file],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT)
    stdout,stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    if output == '':
        print(sv_id + " not found")
        return -1
    else:
        return output.rstrip().split("\n")

def search_gene_info(in_file, gene_id):
    cmd = subprocess.Popen(['zgrep', "-w", gene_id + "\|#Chr" , in_file],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT)
    stdout,stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    if output == '':
        print(gene_id + " not found")
        return -1
    else:
        return output.rstrip().split("\n")

def extract_values(filename):
    cmd = subprocess.Popen(['sort', "-u", filename],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT)
    stdout,stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    if output == '':
        print("no snp was found")
        return -1
    else:
        return output.rstrip()

def standardization(data):
    return data/np.nanstd(data)

def qnorm(beta, p):
    if float(beta) >0:
        return str(-norm.ppf(float(p)/2))
    else:
        return str(norm.ppf(float(p)/2))

final_out = open(egene_file + ".caviar.results", "w")
final_out.write("gene_name\tcaviar_result\n")

for gene_id in egene_list:
    all_cis_markers = search_gene_name(query_file, gene_id)
    sv_list, snp_indel_list = seperate_two_lists(all_cis_markers, ["INS", "DEL"])

    rank_list(sv_list)
    rank_list(snp_indel_list)

    # prepare snp indel SV IDs for plink to extract .bed file for calculating LD
    # prepare zscore based on fastQTL output, by using the qnorm() function
    out = open(folder + gene_id  + ".snp_indel_sv", "w")
    out2= open(folder + gene_id  + ".zscore", "w")
    out.write(sv_list[0][1] + "\n" )
    out2.write(sv_list[0][1] + "\t" + qnorm(sv_list[0][4], sv_list[0][3]) + "\n")
    for i in range(100):
        out.write(snp_indel_list[i][1] + "\n")
        out2.write(snp_indel_list[0][1] + "\t" + qnorm(snp_indel_list[0][4], snp_indel_list[0][3]) + "\n")
    out.close()
    out2.close()

    ## run plink to obtain LD matrix
    cmd = "sh run_plink_for_every_eGene.sh " + folder + gene_id  + ".snp_indel_sv" 
    os.system(cmd)

    # The following is plink.ld
    #  CHR_A         BP_A               SNP_A  CHR_B         BP_B               SNP_B            R 
    #  3    108222323      3_108222323_SV      3    108491648         3_108491648     0.980003 

    ## convert the plink.ld to 2D matrix, as the CAVIAR input
    all_markers = [line.rstrip() for line in open(folder + gene_id  + ".snp_indel_sv").readlines()]
    data = [line.rstrip().split() for line in open(folder + gene_id  + ".snp_indel_sv.plink_format.ld").readlines()[1:]]
    out = open(folder + gene_id  + ".snp_indel_sv.plink_format.ld.2d", "w")
    ld = {}
    for d in data:
        ld[(d[2],d[5])] = d[6]
        ld[(d[5],d[2])] = d[6]
    for m1 in all_markers:   # 2D matrix每一行
        for m2 in all_markers:
            if m1 == m2:
                out.write("1")
            elif ((m1, m2) in ld) :
                out.write(ld[(m1,m2)])    ## 如果有的SNP pair 不在ld 文件中，只能填上0
            else:
                out.write("0")
            if m2 != all_markers[-1]:
                out.write(" ")
            else:
                out.write("\n")
    out.close()

    print("start to run CAVIAR...\n")
    ### run CAVIAR for every gene
    cmd = "sh run_CAVIAR_for_every_eGene.sh " + folder + gene_id  + ".snp_indel_sv " + folder + gene_id  + ".zscore"
    os.system(cmd)

#     #### summarize results in egene_file.results
    if os.path.isfile(folder + gene_id + ".snp_indel_sv.caviar.out_set"):
        variation_id = extract_values(folder + gene_id + ".snp_indel_sv.caviar.out_set") 
        final_out.write(gene_id + "\t" + variation_id + "\n")
        cmd = "rm " + folder + gene_id + ".*"
        os.system(cmd)
    else:
        final_out.write(gene_id + "\tcaviar_failed\n")
        cmd = "rm " + folder + gene_id + ".*"
        os.system(cmd)
final_out.close()
