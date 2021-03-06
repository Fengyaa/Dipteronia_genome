#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd 
import itertools
import gzip
import subprocess
import statistics

###ref https://www.nature.com/articles/ng.3186#Sec8
    

def read_hap_df(vcf_fn, chrom=None, start=None, end=None, samples=None, **kwa): 
    """
    A slightly more advanced vcf parser. Reads in haplotypes from a vcf file. Basically does the same as done in the cells above, but allows the used to specify the range of the genome that should be read in. Also allows to specify which samples should be used.
    Parameters:
    vcf_fn : file path of the VCF to be read
    chrom : specify which chromosome (or scaffold)
    to read from the file
    (only works on bgzipped, tabix-indexed files) default ... read whole file
    start: specify the start nucleotide position (only works if chrom given on bgzipped, tabix-indexed files); default=1
    end: specify the ebd nucleotide position
    (only works if chrom given on bgzipped, tabix-indexed files); default=chrom_end
    samples: list of sample names to read; default ... all samples
    Parameters:
    vcf_fn : file path of the VCF to be read
    chrom : specify which chromosome (or scaffold)
    to read from the file
    (only works on bgzipped, tabix-indexed files) default ... read whole file
    start: specify the start nucleotide position (only works if chrom given on bgzipped, tabix-indexed files); default=1
    end: specify the ebd nucleotide position
    (only works if chrom given on bgzipped, tabix-indexed files); default=chrom_end
    samples: list of sample names to read; default ... all samples
        returns:
        Pandas dataframe of index (chrom, pos)
    and columns (sample, haplotype). are 0 for first and 1 for second
    ###usage example
    # small_hap_df = read_hap_df(vcf_fn, chrom='chr01',start=938478,end=939478)
    # hap_df=read_hap_df(vcf_fn)

    #   hap_df.loc['HIC_ASM_0'].loc[10000:15000]
    """
    # parse header
    header_stream = subprocess.Popen(['bgzip', '-d', vcf_fn, "-c"], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    header_stdin=header_stream.stdout
    # with open(vcf_fn,'rb') as f:
    for line in header_stdin:
        # print(line)
        line = line.decode()
        if line[:6] == '#CHROM':
            vcf_header = line.strip().split('\t') 
            vcf_header[0] = vcf_header[0][1:]
            break
        # determine genomic region to read in
    if chrom is not None:
        assert vcf_fn[-3:] == ".gz", "Only supply chrom if vcf is bgzipped and tabix indexed"
        region = chrom
        if end is not None and start is None:
            start = 0
        if start is not None:
            region += ':' + str(start) 
        if end is not None:
            region += '-' + str(end) 
    else:
        region = None
    # If no specific samples given, use all samples in the VCF
    if samples is None:
        samples = vcf_header[9:]
    # Either use regional input or input whole VCF
    if region is None:
        stdin = vcf_fn 
    else:
        tabix_stream = subprocess.Popen(['tabix', vcf_fn, region], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdin = tabix_stream.stdout
    gen_df = pd.read_csv(stdin, sep='\t',comment='#',names=vcf_header, usecols=['CHROM','POS']+samples, index_col=['CHROM','POS'],dtype=str, **kwa)
    # print(gen_df)
    first_haplotype = gen_df.applymap(lambda s: s.split(":")[0].split('/')[0])
    second_haplotype = gen_df.applymap(lambda s: s.split(":")[0].split('/')[1])
    first_haplotype.columns = pd.MultiIndex.from_product([first_haplotype.columns, [0]])
    second_haplotype.columns = pd.MultiIndex.from_product([second_haplotype.columns, [1]])
    hap_df = pd.concat([first_haplotype, second_haplotype], axis=1).sort_index(axis=1)
    return hap_df
    # return gen_df


popx = [x.strip() for x in open(sys.argv[1], "r")]
popy = [x.strip() for x in open(sys.argv[2], "r")]
site_list=[x.strip() for x in open('../reanalyze1019/Dis_sift_provean.pos',"r")]
vcf_fn = "Dis_all.snpeff_ex_wrong_pol_GT.recode.vcf.gz"

popx_hap_df = read_hap_df(vcf_fn,samples=popx)
popy_hap_df = read_hap_df(vcf_fn,samples=popy)


def rxy_calc(site_list,popx_hap_df,popy_hap_df):
    Lx_not_y=0
    Ly_not_x=0
    for line in site_list:
        x = line.split('\t')
        chr_name = x[0]
        pos = int(x[1])
        gt_x_list = popx_hap_df.loc[(chr_name,pos),:].values.tolist()
        nx = gt_x_list.count("1") + gt_x_list.count("0")
        dx = gt_x_list.count("1")
        gt_y_list = popy_hap_df.loc[(chr_name,pos),:].values.tolist()
        ny = gt_y_list.count("1") + gt_y_list.count("0")
        dy = gt_y_list.count("1")
        # print(dx,nx,dy,ny)
        Lx_not_y += (dx/nx)*(1-dy/ny)
        Ly_not_x += (dy/ny)*(1-dx/nx)

    return(Lx_not_y/Ly_not_x)


window_size=len(site_list)//20
site_lists = []
s=0
for index in range(0,20):
    new_site_list = site_list.copy()
    for i in range(window_size*index,window_size*(index + 1)):
        try:
            new_site_list.pop(i)
        except IndexError:
            continue
    s+=1
    # print(s)
    # print(len(new_site_list))
    site_lists.append(new_site_list)
print(f'total rxy is : {rxy_calc(site_list,popx_hap_df,popy_hap_df)}')
#jackknifing
for s_list in site_lists:
    print(rxy_calc(s_list,popx_hap_df,popy_hap_df))

