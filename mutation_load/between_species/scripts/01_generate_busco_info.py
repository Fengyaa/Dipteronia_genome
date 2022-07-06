#!/usr/bin/env python3

import sys
import pandas as pd 
Did_busco_genes = open('Did_busco_genes.txt')
Did_busco_dict = {}
for line in Did_busco_genes:
    x = line.strip().split(':')
    Did_busco_dict[x[0]] = x[1]

Dis_busco_genes = open('Dis_busco_genes.txt')
Dis_busco_dict = {}
for line in Dis_busco_genes:
    x = line.strip().split(':')
    Dis_busco_dict[x[0]] = x[1]

OG_df = pd.read_csv('SingleCopyOrthologues.tsv',sep = '\t', header = None)

for key in Dis_busco_dict.keys():
    if key in Did_busco_dict.keys():
        dis_gene_name = Dis_busco_dict[key]
        did_gene_name = Did_busco_dict[key]
        if dis_gene_name in OG_df[3].values.tolist() and OG_df[OG_df[3] == dis_gene_name][2].tolist()[0] == did_gene_name:
            OG_name = OG_df[OG_df[3] == dis_gene_name][0].tolist()[0]
            Acy_gene_name = OG_df[OG_df[3] == dis_gene_name][1].tolist()[0]
            print(OG_name,key,did_gene_name,dis_gene_name,Acy_gene_name)
