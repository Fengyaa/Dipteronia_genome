#!/usr/bin/env python3

import sys
import os

def read_fasta(input_fasta):
    align = {} #define a dictionary
    with open(input_fasta, "r") as in_fasta:
        seq = ""
        seqname = ""
        seq_num = 0
        # for x in line:
        for line in in_fasta:
            if line[0] == ">":
                seqname = line.strip(">").strip()
                seq_num = seq_num + 1
                seq = ""
            else:
                seq = seq + line.strip()
            align[seqname] = seq
        return align

Did_pep = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Did.pep")
Did_cds = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Did1.cds")

Dis_pep = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Dis.pep")
Dis_cds = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Dis1.cds")
Did_pep = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Did.pep")
Did_cds = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Did1.cds")

Acy_pep = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Acer.pep")
Acy_cds = read_fasta("/data/xubo/fengyu/dipteronia/01_Genome_Assemblies/Acyan1.cds")


with open('../busco_og_info.txt','r') as in_info:
    for line in in_info:
        x = line.strip().split()
        OG_name = x[0]
        Did_name = x[2]
        Dis_name = x[3]
        Acy_name = x[4]
        with open(f'{OG_name}.cds','w') as out_cds:
            out_cds.write(f'>{Did_name}\n{Did_cds[Did_name]}\n')
            out_cds.write(f'>{Dis_name}\n{Dis_cds[Dis_name]}\n')
            out_cds.write(f'>{Acy_name}\n{Acy_cds[Acy_name]}\n')
        with open(f'{OG_name}.pep','w') as out_pep:
            out_pep.write(f'>{Did_name}\n{Did_pep[Did_name]}\n')
            out_pep.write(f'>{Dis_name}\n{Dis_pep[Dis_name]}\n')
            out_pep.write(f'>{Acy_name}\n{Acy_pep[Acy_name]}\n')    

        os.system(f'clustalw -infile={OG_name}.pep -output=FASTA > /dev/null 2>&1 ')#outputname is ${OG_name}.fasta
        os.system(f'perl /home/xubo/scripts/pal2nal.pl {OG_name}.fasta {OG_name}.cds -output fasta > {OG_name}.codon.fasta')
