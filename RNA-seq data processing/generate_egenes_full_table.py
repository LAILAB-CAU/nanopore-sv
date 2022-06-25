import os
import sys
import subprocess
import numpy as np

egene_file = sys.argv[1]
egene_list= [line.rstrip() for line in open(egene_file).readlines()]

folder = "./"
### for every gene_id, find the 1000 most significant SNPs and small InDels, and the one most signiciant SV
# egene_list= [line.rstrip() for line in open(folder + "/combined.all.chunks.benjamini.txt").readlines()]
gene_expression_file = "DTA73MX-mBsnpREF-200425-expressed_sorted"

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
        d = open(filename, "r").readlines()
        return d[10].rstrip().split("\t")[1], d[1].split("\t")[1], d[1].rstrip().split("\t")[2], d[14].split("\t")[0],d[14].rstrip().split("\t")[1]

def standardization(data):
    return data/np.nanstd(data)

final_out = open(egene_file + ".fastQTL.full_table", "w")
final_out.write("gene_name\tsv_id\tsv_dist\tsv_pvalue\tsv_slope\tsnp_indel_id\tsnp_indel_dist\tsnp_indel_pvalue\tsnp_indel_slope\tvar(FPKM)\tmedian(FPKM)\tmax(FPKM)\n")

for gene_id in egene_list:
    all_cis_markers = search_gene_name(query_file, gene_id)
    sv_list, snp_indel_list = seperate_two_lists(all_cis_markers, ["INS", "DEL"])

    rank_list(sv_list)
    rank_list(snp_indel_list)
    
    print("start to save variant files\n")
    # 输出效应最大的一个SV和一个snp_indel的信息，以及基因的FPKM的信息
    sv_id, sv_dist, sv_pvalue, sv_slope = sv_list[0][1], sv_list[0][2], sv_list[0][3], sv_list[0][4]
    final_out.write(gene_id + "\t" + sv_id + "\t" + sv_dist + "\t" + sv_pvalue + "\t" + sv_slope + "\t")

    snp_indel_id, snp_indel_dist, snp_indel_pvalue, snp_indel_slope = snp_indel_list[0][1], snp_indel_list[0][2], snp_indel_list[0][3], snp_indel_list[0][4]
    final_out.write(snp_indel_id + "\t" + snp_indel_dist + "\t" + snp_indel_pvalue + "\t" + snp_indel_slope + "\t")

    # calculate each egene var(FPKM), median(FPKM), maximum(FPKM)
    gene_info = [item.split("\t") for item in search_gene_info(gene_expression_file, gene_id)]
    express_level = np.array([float(item) for item in gene_info[0][1:]])
    final_out.write(str(np.var(express_level)) + "\t" + str(np.median(express_level)) + "\t" + str(np.max(express_level)) + "\n")
final_out.close()

