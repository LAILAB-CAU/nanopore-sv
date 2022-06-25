import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import os

def normalize_quantiles(M, inplace=False):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")  
    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    if not inplace:
        M = M.copy()
    
    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n
    
    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1
                
        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1
    
    if not inplace:
        return M


def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q
        
        
def get_donors_from_vcf(vcfpath):
    """
    Extract donor IDs from VCF
    """
    with gzip.open(vcfpath) as vcf:
        for line in vcf:
            if line.decode()[:2]=='##': continue
            break
    return line.decode().strip().split('\t')[9:]

def get_donors_from_fpkm(fpkmpath):
    """
    Extract donor IDs from FPKM table
    """
    return open(fpkmpath).readline().rstrip().split(',')[1:]


def normalize_expression(expression_df, expression_threshold=0.1,min_samples=10):
    """
    Genes are thresholded based on the following expression rules:
      >=min_samples with >expression_threshold expression values
    """
    donor_ids = ['-'.join(i.split('-')[:2]) for i in expression_df.columns]
    
    # expression thresholds
    mask = (np.sum(expression_df>expression_threshold,axis=1)>=min_samples).values
    
    # apply normalization
    M = normalize_quantiles(expression_df.loc[mask].values, inplace=False)
    R = inverse_quantile_normalization(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index)    
    # quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index)
    return quant_std_df #, quant_df

    
def read_gct(gct_file, donor_ids):
    """
    Load GCT as DataFrame
    """    
    df = pd.read_csv(gct_file, sep=',', skiprows=0, index_col=0)
    # df.drop('gene', axis=1, inplace=True) 
    df.index.name = 'gene_id'
    return df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('expression_gct', help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('--expression_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least min_samples')
    parser.add_argument('--min_samples', type=np.int32, default=10, help='Minimum number of samples that must satisfy thresholds')

    args = parser.parse_args()
    
    print('Generating normalized expression files ... ', end='', flush=True)
    donor_ids = get_donors_from_fpkm(args.expression_gct)
    expression_df = read_gct(args.expression_gct, donor_ids)

    quant_std_df = normalize_expression(expression_df, expression_threshold=args.expression_threshold, min_samples=args.min_samples)

    # for consistency with v6/v6p pipeline results, write unsorted expression file for PEER factor calculation
    quant_std_df.to_csv(args.prefix+'.expression.txt', sep='\t')