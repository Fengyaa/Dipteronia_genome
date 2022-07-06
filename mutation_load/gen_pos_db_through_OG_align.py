#!/usr/bin/env python3
import os
from Bio import SeqIO
import pandas as pd
import re
import sys
from multiprocessing import Pool

Dis_cds_gff = pd.read_table("../Dis.ogcds.gff", header=None,sep="\t")
Did_cds_gff = pd.read_table("../Did.ogcds.gff", header=None,sep="\t")
Dis_cds_gff["gene_ID"]= Dis_cds_gff.apply(lambda x : re.findall(r"ID=cds\.(.+?)\;", x[8])[0],axis=1)
Did_cds_gff["gene_ID"]= Did_cds_gff.apply(lambda x : re.findall(r"ID=cds\.(.+?)\;", x[8])[0],axis=1)
Dis_genome_fasta = SeqIO.index("/Users/yufeng/Desktop/Dipteronia_genomes/Dipteornia_sinensis/00.Genome/genome.fa","fasta")
Did_genome_fasta = SeqIO.index("/Users/yufeng/Desktop/Dipteronia_genomes/Dipteronia_dyeriana/00.Genome/genome.fa","fasta")

base_tras = {
    "A":"T",
    "G":"C",
    "C":"G",
    "T":"A"}

def get_gene_cds_range_list(gene_name,cds_gff):
    gene_gff = cds_gff[cds_gff["gene_ID"] == gene_name]
    cds_directions = gene_gff[6].tolist()
    # print(gene_gff[0].tolist())
    chrom = gene_gff[0].tolist()[0]
    assert len(set(cds_directions)) == 1 , "invalid cds gff infomations!"
    cds_direction = cds_directions[0]
    base_pos_list = []
    for index, row in gene_gff.iterrows():
        #real position in row3 vs. index position starts from 0
        i_cds_list = list(range(row[3]-1,row[4]))
        base_pos_list += i_cds_list
    # if cds_direction == "-":
    #     base_pos_list.reverse()
    base_pos_list.sort()
    # print(base_pos_list)
    return [chrom,cds_direction,base_pos_list]


def extract_cds(d_base_pos_list,genome_fasta_dict):
    chrom = d_base_pos_list[0]
    direction = d_base_pos_list[1]
    position_list = d_base_pos_list[2]
    chrom_seq = str(genome_fasta_dict[chrom].seq)
    base_list = []
    if direction == "+":
        for i in position_list:
            base_list.append(chrom_seq[i])
    else:
        position_list.reverse()
        for i in position_list:
            base_list.append(base_tras[chrom_seq[i]])
    cds_seq = "".join(base_list)
    return cds_seq

def join_positions(OG_name):
    cds_align = SeqIO.to_dict(SeqIO.parse(OG_name + '.codon.fasta', "fasta"))
    Did_gene_name = list(cds_align.keys())[0]
    Did_gene_seq = str(cds_align[Did_gene_name].seq)
    Dis_gene_name = list(cds_align.keys())[1]
    Dis_gene_seq = str(cds_align[Dis_gene_name].seq)
    Dis_cds_pos = get_gene_cds_range_list(Dis_gene_name.replace('Dis','evm.model'),Dis_cds_gff)
    Did_cds_pos = get_gene_cds_range_list(Did_gene_name.replace('Did','evm.model'),Did_cds_gff)
    # print(extract_cds(Dis_cds_pos,Dis_genome_fasta)[0:-3] == Dis_gene_seq.replace('-',''))
    # print(Dis_gene_seq.replace('-',''))
    # print(extract_cds(Did_cds_pos,Did_genome_fasta)[0:-3])
    # print(Did_gene_seq.replace('-',''))
    assert extract_cds(Dis_cds_pos,Dis_genome_fasta)[0:-3] == Dis_gene_seq.replace('-',''), 'Dis extracted seq does not equal coding sequences'
    assert extract_cds(Did_cds_pos,Did_genome_fasta)[0:-3] == Did_gene_seq.replace('-',''), 'Did extracted seq does not equal coding sequences'
    gene_pos_df = pd.DataFrame(columns=['align_index','Dis_chrom','Dis_chr_pos','Dis_gene_name','Dis_cds_basepos','Dis_ref','Did_chrom','Did_chr_pos','Did_gene_name','Did_cds_basepos','Did_ref','Orthogroup'])
    for i in range(len(Did_gene_seq)):
        Dis_cds_basepos = len(Dis_gene_seq[0:i]) - Dis_gene_seq[0:i].count('-')
        Did_cds_basepos = len(Did_gene_seq[0:i]) - Did_gene_seq[0:i].count('-')
        pos_inf = {
        'align_index':i,
        'Dis_chrom':Dis_cds_pos[0],
        'Dis_chr_pos':Dis_cds_pos[2][Dis_cds_basepos],
        'Dis_gene_name':Dis_gene_name,
        'Dis_cds_basepos':Dis_cds_basepos,
        'Dis_ref':Dis_gene_seq[i],
        'Did_chrom':Did_cds_pos[0],
        'Did_chr_pos':Did_cds_pos[2][Did_cds_basepos],
        'Did_gene_name':Did_gene_name,
        'Did_cds_basepos':Did_cds_basepos,
        'Did_ref':Did_gene_seq[i],
        'Orthogroup': OG_name,
        }
        if pos_inf['Did_ref'] == '-':
            pos_inf['Did_chrom'] = '-'
            pos_inf['Did_chr_pos'] = '-'
            pos_inf['Did_cds_basepos'] = '-'
        if pos_inf['Dis_ref'] == '-':
            pos_inf['Dis_chrom'] = '-'
            pos_inf['Dis_chr_pos'] = '-'
            pos_inf['Dis_cds_basepos'] = '-'
        gene_pos_df = gene_pos_df.append(pos_inf,ignore_index=True)

    return gene_pos_df

def write_result(OG):
    if os.path.exists(f'{OG}.pos_db.csv') == False:
        try:
            joined_df = join_positions(OG)
            joined_df.to_csv(f'{OG}.pos_db.csv',index=False)
        except:
            next

# testdf = join_positions('OG0015648')
# testdf.to_csv(f'../test_df',index=False)

OGs = [x.strip().split()[0] for x in open('../SingleCopyOrthologues.tsv','r')]
for OG in OGs:
    joined_df = join_positions(OG)
    joined_df.to_csv(f'{OG}.pos_db.csv',index=False)

# def run_complex_operations(operation, input, pool):
#     pool.map(operation, input)

# processes_count = 20

# if __name__ == '__main__':
#     processes_pool = Pool(processes_count)
#     run_complex_operations(write_result, OGs, processes_pool)   
