#!/usr/bin/env python3

#extract cds from vcf to fasta
#note: positions generated here are started from zero
#Input files: vcf_file, cdsgfffile, genomefasta
import sys
import pandas as pd 
import numpy as np
import itertools
import gzip
import subprocess
import re
from Bio import SeqIO
import os 

base_tras = {
    "A":"T",
    "G":"C",
    "C":"G",
    "T":"A",
    "-":"-"}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

NUCLEOTIDE_BASE = {
    "DNA": ["A", "T", "C", "G"],
    "RNA": ["A", "U", "C", "G"]
}

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

def create_gt_dict(sample_list,gt):
    d = {}
    for i in sample_list:
        d[i] = gt
    return d

def trans_gt_dict(d):
    trans_dict = {'0/0':'1/1','1/1':'0/0','0/1':'1/0','1/0':'0/1','.':'.','./.':'./.'}
    for i, j in d.items():
        d[i] = trans_dict[j]
    return(d)


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
    gen_df = pd.read_csv(stdin,sep='\t',comment='#',header=None,names=vcf_header, index_col=['CHROM','POS'], **kwa)
    ref_alt_df = gen_df.loc[:, ['REF', 'ALT']]
    geno_df = gen_df.iloc[:, 7:].applymap(lambda s: s.split(":")[0])
    hap_df = pd.concat([ref_alt_df,geno_df],axis=1)
    return hap_df
# print(read_hap_df('Did.ogcds_rmdp_sorted.vcf.gz'))

def join_positions(OG_name,Dis_hap_df,Did_hap_df):
    cds_align = SeqIO.to_dict(SeqIO.parse(OG_name + '.codon.fasta', "fasta"))
    ancestral_seq = open(OG_name + '.anc.cds').readlines()[0]
    Dis_samples = Dis_hap_df.columns.tolist()[2:]
    Did_samples = Did_hap_df.columns.tolist()[2:]
    #get gene names for Did and Dis
    for x in list(cds_align.keys()):
        if "Dis" in x:
            Dis_gene_name = x
        elif 'Did' in x:
            Did_gene_name = x
        else:
            continue
    #get sequences 
    Did_gene_seq = str(cds_align[Did_gene_name].seq)
    Dis_gene_seq = str(cds_align[Dis_gene_name].seq)

    #get cds postision dict
    Dis_cds_pos = get_gene_cds_range_list(Dis_gene_name.replace('Dis','evm.model'),Dis_cds_gff)
    Did_cds_pos = get_gene_cds_range_list(Did_gene_name.replace('Did','evm.model'),Did_cds_gff)
    # print(extract_cds(Dis_cds_pos,Dis_genome_fasta)[0:-3] == Dis_gene_seq.replace('-',''))
    # print(Dis_gene_seq.replace('-',''))
    # print(extract_cds(Did_cds_pos,Did_genome_fasta)[0:-3])
    # print(Did_gene_seq.replace('-',''))
    
    assert extract_cds(Dis_cds_pos,Dis_genome_fasta)[0:-3] == Dis_gene_seq.replace('-',''), 'Dis extracted seq does not equal coding sequences'
    assert extract_cds(Did_cds_pos,Did_genome_fasta)[0:-3] == Did_gene_seq.replace('-',''), 'Did extracted seq does not equal coding sequences'
    # gene_pos_df = pd.DataFrame(columns=['align_index','Dis_chrom','Dis_chr_pos','Dis_gene_name','Dis_gene_direction','Dis_cds_basepos','Dis_ref','Did_chrom','Did_chr_pos','Did_gene_name','Did_gene_direction','Did_cds_basepos','Did_ref','Orthogroup'])
    gene_vcf_df = pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + Dis_samples + Did_samples)

    for i in range(len(Did_gene_seq)):
        Dis_cds_basepos = len(Dis_gene_seq[0:i]) - Dis_gene_seq[0:i].count('-') ## start from zero
        Did_cds_basepos = len(Did_gene_seq[0:i]) - Did_gene_seq[0:i].count('-')

        ancestral_base = ancestral_seq[i]
        if ancestral_base == '-':
            continue

        pos_inf = {
        'align_index':i,
        'Dis_chrom':Dis_cds_pos[0],
        'Dis_chr_pos':Dis_cds_pos[2][Dis_cds_basepos],
        'Dis_gene_name':Dis_gene_name,
        'Dis_gene_direction':Dis_cds_pos[1],
        'Dis_cds_basepos':Dis_cds_basepos,
        'Dis_ref':Dis_gene_seq[i],
        'Did_chrom':Did_cds_pos[0],
        'Did_chr_pos':Did_cds_pos[2][Did_cds_basepos],
        'Did_gene_name':Did_gene_name,
        'Did_gene_direction':Did_cds_pos[1],
        'Did_cds_basepos':Did_cds_basepos,
        'Did_ref':Did_gene_seq[i],
        'Orthogroup': OG_name,
        }

        if pos_inf['Did_ref'] == '-':
            pos_inf['Did_chrom'] = '-'
            pos_inf['Did_chr_pos'] = '-'
            pos_inf['Did_cds_basepos'] = '-'
            Did_ref = '.'
            Did_alt = '.'
            did_gt_dict = create_gt_dict(Did_samples,'.')
        else:
            Did_chr_pos = int(pos_inf['Did_chr_pos']) + 1
            if (pos_inf['Did_chrom'],Did_chr_pos) in Did_hap_df.index:
                cur_did_hap_df = Did_hap_df.loc[(pos_inf['Did_chrom'],Did_chr_pos)]
                did_gt_dict = cur_did_hap_df.to_dict()
                Did_ref = did_gt_dict.pop('REF')
                Did_alt = did_gt_dict.pop('ALT')
                if pos_inf['Did_gene_direction'] == '-':
                    Did_ref = base_tras[Did_ref]
                    Did_alt = base_tras[Did_alt]
                    assert pos_inf['Did_ref'] == Did_ref,f'Did ref does not consistent in {pos_inf["Did_chrom"]} {Did_chr_pos} {OG_name} {i}'
            else:
                Did_ref = pos_inf['Did_ref']
                Did_alt = '.'
                did_gt_dict = create_gt_dict(Did_samples,'0/0')
            
        if pos_inf['Dis_ref'] == '-':
            pos_inf['Dis_chrom'] = '-'
            pos_inf['Dis_chr_pos'] = '-'
            pos_inf['Dis_cds_basepos'] = '-'
            Dis_ref = '.'
            Dis_alt = '.'
            dis_gt_dict = create_gt_dict(Dis_samples,'.')
        else:
            Dis_chr_pos = int(pos_inf['Dis_chr_pos']) + 1
            if (pos_inf['Dis_chrom'],Dis_chr_pos) in Dis_hap_df.index:
                cur_dis_hap_df = Dis_hap_df.loc[(pos_inf['Dis_chrom'],Dis_chr_pos)]
                dis_gt_dict = cur_dis_hap_df.to_dict()
                Dis_ref = dis_gt_dict.pop('REF')
                Dis_alt = dis_gt_dict.pop('ALT')
                if pos_inf['Dis_gene_direction'] == '-':
                    Dis_ref = base_tras[Dis_ref]
                    Dis_alt = base_tras[Dis_alt]
                    assert pos_inf['Dis_ref'] == Dis_ref,f'Dis ref does not consistent in {pos_inf["Dis_chrom"]} {Dis_chr_pos} {OG_name} {i}'
            else:
                Dis_ref = pos_inf['Dis_ref']
                Dis_alt = '.'
                dis_gt_dict = create_gt_dict(Dis_samples,'0/0')

        base_set_list = list(set([ancestral_base,Did_ref,Did_alt,Dis_ref,Dis_alt]))
        if '.' in base_set_list:
            base_set_list.remove('.')
        if '-' in base_set_list:
            base_set_list.remove('-')

        if len(base_set_list) != 2:
            continue
        else:
            base_set_list.remove(ancestral_base)
            alt_base = base_set_list[0]

            if Did_ref == alt_base:
                did_gt_dict = trans_gt_dict(did_gt_dict)
            if Dis_ref == alt_base:
                dis_gt_dict = trans_gt_dict(dis_gt_dict)

            #get ancestral codon and peptide
            ancestral_codon = 'none'
            ancestral_codon = ancestral_seq[ (i - i % 3):(i - i % 3 + 3) ]
            ancestral_pep = DNA_Codons[ancestral_codon]
            pep_pos = i // 3 + 1
            base_codon_pos = i % 3
            alt_codon = ancestral_codon[0:base_codon_pos] + alt_base + ancestral_codon[base_codon_pos+1:]
            alt_pep = DNA_Codons[alt_codon]
            # print(OG_name,i,ancestral_codon,ancestral_pep,pep_pos,base_codon_pos,alt_base,alt_codon,alt_pep)
            if alt_pep == '_':
                Effect = "LOF"
            if ancestral_pep == alt_pep :
                Effect = 'synonymous'
            else:
                Effect = ancestral_pep + str(pep_pos) + alt_pep
            # print(Effect)
            pos_vcf_info = {
            '#CHROM': OG_name,
            'POS': i+1,
            'ID': '.',
            'REF': ancestral_base,
            'ALT': alt_base,
            'QUAL':'.',
            'FILTER':'.',
            'INFO':Effect,
            'FORMAT':'GT',
            }
            pos_vcf_info.update(did_gt_dict)
            pos_vcf_info.update(dis_gt_dict)

        gene_vcf_df = gene_vcf_df.append(pos_vcf_info,ignore_index=True)
    return gene_vcf_df

def write_result(OG):
    if os.path.exists(f'{OG}.busco.vcf') == False:
        try:
            joined_df = join_positions(OG,Dis_hap_df,Did_hap_df)
            joined_df.to_csv(f'{OG}.busco.vcf',sep="\t",index=False)
        except:
            print(f'{OG} failed')
            next

Dis_cds_gff = pd.read_table("../Dis.ogcds.gff", header=None,sep="\t")
Did_cds_gff = pd.read_table("../Did.ogcds.gff", header=None,sep="\t")
Dis_cds_gff["gene_ID"]= Dis_cds_gff.apply(lambda x : re.findall(r"ID=cds\.(.+?)\;", x[8])[0],axis=1)
Did_cds_gff["gene_ID"]= Did_cds_gff.apply(lambda x : re.findall(r"ID=cds\.(.+?)\;", x[8])[0],axis=1)
Dis_genome_fasta = SeqIO.index("/Users/yufeng/Desktop/Dipteronia_genomes/Dipteornia_sinensis/00.Genome/genome.fa","fasta")
Did_genome_fasta = SeqIO.index("/Users/yufeng/Desktop/Dipteronia_genomes/Dipteronia_dyeriana/00.Genome/genome.fa","fasta")
Did_hap_df = read_hap_df('../Did.ogcds_rmdp_sorted.vcf.gz')
Dis_hap_df = read_hap_df('../Dis.ogcds_rmdp_sorted.vcf.gz')


# res = join_positions('OG0016949',Dis_hap_df,Did_hap_df)
# write_result('OG0014872')
# print(res)
# testdf = join_positions('OG0015648')
# testdf.to_csv(f'test_df',index=False)
for line in open('OG_names.txt'):
    OG_name = line.strip()
    write_result(OG_name)
