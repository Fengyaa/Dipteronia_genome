#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

Dis_gene_seqs = SeqIO.to_dict(SeqIO.parse('Dis_busco_scngs/busco_scngs.fa', "fasta"))
Did_gene_seqs = SeqIO.to_dict(SeqIO.parse('Did_busco_scngs/busco_scngs.fa', "fasta"))

def read_sift_score(gene_name,species_name,gene_seq_dict):
    x = pd.read_csv(f'{species_name}_busco_scngs/{gene_name}.SIFTprediction',sep = '\s+', skiprows=4)
    gene_fasta =  str(gene_seq_dict[gene_name].seq)
    # assert len(gene_fasta) == len(x), f'{gene_name},gene_length = {len(gene_fasta)}, db_length = {len(x)}'
    i = 0
    num_del = 0
    for index,row in x.iterrows():
        y = row[gene_fasta[i]]
        i += 1
        if i > len(gene_fasta) - 1:
            break
        if float(y) < 0.05 :
            num_del += 1
    return num_del
# read_sift_score('evm.model.HIC_ASM_9.295','Dis',Dis_gene_seqs)

print('OG_name busco_name Dis_gene_name Dis_del_num Did_gene_name Did_del_num')
for line in open('busco_og_info.txt'):
    x = line.strip().split(' ')
    OG_name = x[0]
    busco_name = x[1]
    Did_gene_name = x[2].replace('Did','evm.model')
    Dis_gene_name = x[3].replace('Dis','evm.model')
    Dis_del_num = read_sift_score(Dis_gene_name,'Dis',Dis_gene_seqs)
    Did_del_num = read_sift_score(Did_gene_name,'Did',Did_gene_seqs)
    print(OG_name,busco_name,Dis_gene_name,Dis_del_num,Did_gene_name,Did_del_num)
