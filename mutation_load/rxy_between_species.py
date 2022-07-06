#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd 
import itertools
import gzip
import subprocess


def read_hap_df(vcf_fn, chrom=None, start=None, end=None, samples=None, **kwa): 
    """
    ###usage example
    # small_hap_df = read_hap_df(vcf_fn, chrom='chr01',start=938478,end=939478)
    # hap_df=read_hap_df(vcf_fn)
    #   hap_df.loc['HIC_ASM_0'].loc[10000:15000]
    """
    # parse header
    header_stream = subprocess.Popen(['bgzip', '-d', vcf_fn, "-c"], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    header_stdin=header_stream.stdout
    with open(vcf_fn,'r') as f:
        # for line in f:
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
    # gen_df = pd.read_csv(stdin, sep='\t',comment='#',names=vcf_header, usecols=['CHROM','POS','REF','ALT']+samples, index_col=['CHROM','POS'], **kwa)
    gen_df = pd.read_csv(stdin,sep='\t',comment='#',header=None,names=vcf_header,    index_col=['CHROM','POS'], **kwa)
    ref_alt_df = gen_df.loc[:, ['REF', 'ALT']]
    # first_haplotype = gen_df.iloc[:, 7:].applymap(lambda s: int(s.split(":")[0].split('/')[0])) 
    # second_haplotype = gen_df.iloc[:, 7:].applymap(lambda s: int(s.split(":")[0].split('/')[1]))
    first_haplotype = gen_df.iloc[:, 7:].applymap(lambda s: s.split(":")[0].split('/')[0])
    second_haplotype = gen_df.iloc[:, 7:].applymap(lambda s: s.split(":")[0].split('/')[1])
    first_haplotype.columns = pd.MultiIndex.from_product([first_haplotype.columns, [0]])
    second_haplotype.columns = pd.MultiIndex.from_product([second_haplotype.columns, [1]])
    # hap_df = pd.concat([first_haplotype, second_haplotype], axis=1).sort_index(axis=1)
    hap_df = pd.concat([ref_alt_df,first_haplotype,second_haplotype], axis=1)
    return hap_df

base_tras = {
    "A":"T",
    "G":"C",
    "C":"G",
    "T":"A",
    "-":"-"}

Did_hap_df = read_hap_df('Did.ogcds_rmdp_sorted.vcf.gz')
Dis_hap_df = read_hap_df('Dis.ogcds_rmdp_sorted.vcf.gz')
position_db = pd.read_csv('OG_dbs.csv',index_col = ['Dis_chrom','Dis_chr_pos']).sort_index()
# position_db = pd.read_csv('OG_dbs.csv',index_col = ['Did_chrom','Did_chr_pos']).sort_index()
# Did_site_list = 'Did_sift_provean.pos' #del
# Dis_site_list = 'Dis_sift_provean.pos'#del
Did_site_list = 'Did_all.ex_wrong_pol.syn.pos'#syn
Dis_site_list = 'Dis_all.ex_wrong_pol.syn.pos'#syn
Did_cds_direction = pd.read_csv("Did.ogcds.direction.txt",sep = ' ',header = None).set_index(0).to_dict()[1]
Dis_cds_direction = pd.read_csv("Dis.ogcds.direction.txt",sep = ' ',header = None).set_index(0).to_dict()[1]

join_del_site = pd.DataFrame(columns=['Dis_chrom','Dis_chr_pos','Dis_gene_name','Dis_ref','Dis_ref_count','Dis_alt','Dis_alt_count','Did_chrom','Did_chr_pos','Did_gene_name','Did_ref','Did_ref_count','Did_alt','Did_alt_count','gene_direction','ref_genome'])


for line in open(Dis_site_list,"r"): #use deleterious sites of dis as query
    # x = line.strip().split("\t")
    x = line.strip().split(" ")
    Dis_chrom = x[0]
    Dis_chr_pos = x[1] #position in the vcf file,1 based
    if (Dis_chrom,str(int(Dis_chr_pos)-1)) in position_db.index:
        cur_pos_db = position_db.loc[(Dis_chrom,str(int(Dis_chr_pos)-1))] ##position in position_db is 0 based
        Dis_gene_name =  cur_pos_db['Dis_gene_name'].values.tolist()[0]
        Dis_gene_direction = Dis_cds_direction[Dis_gene_name]
        Did_gene_name =  cur_pos_db['Did_gene_name'].values.tolist()[0]
        Did_gene_direction = Did_cds_direction[Did_gene_name]
        cur_dis_hap_df = Dis_hap_df.loc[(Dis_chrom,int(Dis_chr_pos))]
        Dis_ref = cur_dis_hap_df['REF']
        Dis_alt = cur_dis_hap_df['ALT']
        Dis_allele_list = cur_dis_hap_df.values.tolist()[2:]
        Dis_ref_count = Dis_allele_list.count("0")
        Dis_alt_count = Dis_allele_list.count("1")
        ##extract information of Did
        Did_chrom = cur_pos_db['Did_chrom'].values.tolist()[0]
        Did_db_pos = cur_pos_db['Did_chr_pos'].values.tolist()[0]
        if Did_db_pos != '-':
            Did_chr_pos = int(Did_db_pos) + 1
            if (Did_chrom,Did_chr_pos) in Did_hap_df.index:
                cur_did_hap_df = Did_hap_df.loc[(Did_chrom,Did_chr_pos)]
                Did_ref = cur_did_hap_df['REF']
                Did_alt = cur_did_hap_df['ALT']
                Did_allele_list = cur_did_hap_df.values.tolist()[2:]
                Did_ref_count = Did_allele_list.count("0")
                Did_alt_count = Did_allele_list.count("1")
            else: 
                #treat this locus in did as invariant site
                if Did_gene_direction == "+" :
                    Did_ref = cur_pos_db['Did_ref'].values.tolist()[0]
                else:
                    Did_ref = base_tras[cur_pos_db['Did_ref'].values.tolist()[0]]
                Did_alt = '-'
                Did_ref_count = 50
                Did_alt_count = 0
            cur_pos_information = {
            'Dis_chrom':Dis_chrom,
            'Dis_chr_pos':Dis_chr_pos,
            'Dis_gene_name':Dis_gene_name,
            'Dis_gene_direction':Dis_gene_direction,
            'Dis_ref':Dis_ref,
            'Dis_ref_count':Dis_ref_count,
            'Dis_alt':Dis_alt,
            'Dis_alt_count':Dis_alt_count,
            'Did_chrom':Did_chrom,
            'Did_chr_pos':Did_chr_pos,
            'Did_gene_name':Did_gene_name,
            'Did_gene_direction':Did_gene_direction,
            'Did_ref':Did_ref,
            'Did_ref_count':Did_ref_count,
            'Did_alt':Did_alt,
            'Did_alt_count':Did_alt_count,
            'ref_genome':'Dis'
            }
            join_del_site = join_del_site.append(cur_pos_information,ignore_index=True)
join_del_site = join_del_site.drop_duplicates()
join_del_site.to_csv('Dis_syn_pos_db.csv',index=False)

position_db = pd.read_csv('OG_dbs.csv',index_col = ['Did_chrom','Did_chr_pos']).sort_index()
for line in open(Did_site_list,"r"): #use deleterious sites of did as query
    x = line.strip().split(' ')
    Did_chrom = x[0]
    Did_chr_pos = x[1] #position in the vcf file,1 based
    if (Did_chrom,str(int(Did_chr_pos)-1)) in position_db.index:
        cur_pos_db = position_db.loc[(Did_chrom,str(int(Did_chr_pos)-1))] ##position in position_db is 0 based
        Dis_gene_name =  cur_pos_db['Dis_gene_name'].values.tolist()[0]
        Dis_gene_direction = Dis_cds_direction[Dis_gene_name]
        Did_gene_name =  cur_pos_db['Did_gene_name'].values.tolist()[0]
        Did_gene_direction = Did_cds_direction[Did_gene_name]
        cur_did_hap_df = Did_hap_df.loc[(Did_chrom,int(Did_chr_pos))]
        Did_ref = cur_did_hap_df['REF']
        Did_alt = cur_did_hap_df['ALT']
        Did_allele_list = cur_did_hap_df.values.tolist()[2:]
        Did_ref_count = Did_allele_list.count("0")
        Did_alt_count = Did_allele_list.count("1")
        ##extract information of Dis
        Dis_chrom = cur_pos_db['Dis_chrom'].values.tolist()[0] 
        Dis_db_pos = cur_pos_db['Dis_chr_pos'].values.tolist()[0]
        if Dis_db_pos != '-':
            Dis_chr_pos = int(Dis_db_pos) + 1
            if (Dis_chrom,Dis_chr_pos) in Dis_hap_df.index:
                cur_dis_hap_df = Dis_hap_df.loc[(Dis_chrom,Dis_chr_pos)]
                Dis_ref = cur_dis_hap_df['REF']
                Dis_alt = cur_dis_hap_df['ALT']
                Dis_allele_list = cur_dis_hap_df.values.tolist()[2:]
                Dis_ref_count = Dis_allele_list.count("0")
                Dis_alt_count = Dis_allele_list.count("1")
            else: 
                #treat this locus in did as invariant site
                if Dis_gene_direction == "+" :
                    Dis_ref = cur_pos_db['Dis_ref'].values.tolist()[0]
                else:
                    Dis_ref = base_tras[cur_pos_db['Dis_ref'].values.tolist()[0]]
                Dis_alt = '-'
                Dis_ref_count = 108
                Dis_alt_count = 0
            cur_pos_information = {
            'Dis_chrom':Dis_chrom,
            'Dis_chr_pos':Dis_chr_pos,
            'Dis_gene_name':Dis_gene_name,
            'Dis_gene_direction':Dis_gene_direction,
            'Dis_ref':Dis_ref,
            'Dis_ref_count':Dis_ref_count,
            'Dis_alt':Dis_alt,
            'Dis_alt_count':Dis_alt_count,
            'Did_chrom':Did_chrom,
            'Did_chr_pos':Did_chr_pos,
            'Did_gene_name':Did_gene_name,
            'Did_gene_direction':Did_gene_direction,
            'Did_ref':Did_ref,
            'Did_ref_count':Did_ref_count,
            'Did_alt':Did_alt,
            'Did_alt_count':Did_alt_count,
            'ref_genome':'Did'
            }
            join_del_site = join_del_site.append(cur_pos_information,ignore_index=True)

join_del_site = join_del_site.drop_duplicates()
join_del_site.to_csv('Did_syn_pos_db.csv',index=False)


def calc_rxy():
    base_tras = {
        "A":"T",
        "G":"C",
        "C":"G",
        "T":"A",
        "-":"-"}

    Did_join_del_site = pd.read_csv('Did_syn_pos_db.csv').sample(n=2000)
    Dis_join_del_site = pd.read_csv('Dis_syn_pos_db.csv').sample(n=2000)
    Lx_not_y = 0
    Ly_not_x = 0
    num_of_site_used = 0
    for index, row in Did_join_del_site.iterrows():
        Dis_ref = row['Dis_ref']
        Dis_alt = row['Dis_alt'] 
        Did_ref = row['Did_ref']
        Did_alt = row['Did_alt']
        # print(Dis_ref,Dis_alt,Did_ref,Did_alt,row['Dis_gene_direction'],row['Did_gene_direction'])
        if row['Dis_gene_direction'] == '-':
            Dis_ref = base_tras[Dis_ref]
            Dis_alt = base_tras[Dis_alt]
        if row['Did_gene_direction'] == '-':
            Did_ref = base_tras[Did_ref]
            Did_alt = base_tras[Did_alt]
        # print([Dis_ref,Dis_alt,Did_ref,Did_alt])
        #check if there's only two allele in the locus
        allele_list = [Dis_ref,Dis_alt,Did_ref,Did_alt]
        if '-' in allele_list:
            allele_list.remove('-')
        # print(allele_list,[row['Dis_ref_count'],row['Dis_alt_count'],row['Did_ref_count'],row['Did_alt_count']])
        if len(set(allele_list)) <= 2 :
            if Did_ref == Dis_ref:#this includes the case that Did_alt == '-' or Dis_alt == '-'
                ## calculate rxy for small vs. large pop, x is did, y is dis
                dx = row['Did_alt_count'] 
                nx = dx + row['Did_ref_count']
                dy = row['Dis_alt_count']
                ny = dy + row['Dis_ref_count']
            elif Did_alt == Dis_ref or Did_ref == Dis_alt:
                dx = row['Did_ref_count'] 
                nx = dx + row['Did_alt_count']
                dy = row['Dis_alt_count']
                ny = dy + row['Dis_ref_count']
            else:
                continue
            Lx_not_y += (dx/nx)*(1-dy/ny)
            Ly_not_x += (dy/ny)*(1-dx/nx)
            num_of_site_used += 1
            # if num_of_site_used > 1000:
            #     break

    for index, row in Dis_join_del_site.iterrows():
        Dis_ref = row['Dis_ref']
        Dis_alt = row['Dis_alt'] 
        Did_ref = row['Did_ref']
        Did_alt = row['Did_alt']
        # print(Dis_ref,Dis_alt,Did_ref,Did_alt,row['Dis_gene_direction'],row['Did_gene_direction'])
        if row['Dis_gene_direction'] == '-':
            Dis_ref = base_tras[Dis_ref]
            Dis_alt = base_tras[Dis_alt]
        if row['Did_gene_direction'] == '-':
            Did_ref = base_tras[Did_ref]
            Did_alt = base_tras[Did_alt]
        # print([Dis_ref,Dis_alt,Did_ref,Did_alt])
        #check if there's only two allele in the locus
        allele_list = [Dis_ref,Dis_alt,Did_ref,Did_alt]
        if '-' in allele_list:
            allele_list.remove('-')
        # print(allele_list,[row['Dis_ref_count'],row['Dis_alt_count'],row['Did_ref_count'],row['Did_alt_count']])
        if len(set(allele_list)) <= 2 :
            if Did_ref == Dis_ref:#this includes the case that Did_alt == '-' or Dis_alt == '-'
                ## calculate rxy for small vs. large pop, x is did, y is dis
                dx = row['Did_alt_count'] 
                nx = dx + row['Did_ref_count']
                dy = row['Dis_alt_count']
                ny = dy + row['Dis_ref_count']
            elif Did_alt == Dis_ref or Did_ref == Dis_alt:
                dx = row['Did_ref_count'] 
                nx = dx + row['Did_alt_count']
                dy = row['Dis_alt_count']
                ny = dy + row['Dis_ref_count']
            else:
                continue
            Lx_not_y += (dx/nx)*(1-dy/ny)
            Ly_not_x += (dy/ny)*(1-dx/nx)
            num_of_site_used += 1
            # if num_of_site_used > 2000:
            #     break

    print(Lx_not_y/Ly_not_x)
    # print(f'a total of {num_of_site_used} used')
for i in range(0,20):
    calc_rxy()
# print(f'a total of {num_of_site_used} used')
# print(f'rxy value is {Lx_not_y/Ly_not_x}')

