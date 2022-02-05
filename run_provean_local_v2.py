#!/usr/bin/env python3
# Usage: python run_provean.py pep_file SIFTannotations.xls

# We suggest using a cutoff of -2.5 for the PROVEAN score when using the
# NCBI nr protein database released in August 2011. That is, consider a score
# higher than -2.5 to be neutral (tolerated) and that lower than or equal to
# -2.5 to be deleterious (damaging). The PROVEAN scores and optimal cutoff
# may slightly vary with different versions of nr database because the scores
# are computed based on the homologs in the DB. More detailed information on
# PROVEAN scores can be found at http://provean.jcvi.org/about.php

import sys
import os
import pandas as pd
from Bio import SeqIO,AlignIO
import re
from datetime import datetime
import logging
import traceback

input_pep = sys.argv[1]
input_vep = sys.argv[2]
sp_name = input_pep.split(".")[0]
provean_path = '/mnt/sdc1/luruisen/fengyu/provean/bin/provean.sh'
# finished_set= set([x.strip().split("\t")[0] for x in open("Dis_missense_finished_genes.txt","r")])

logging.basicConfig(filename= sp_name + "provean_log.txt",
                    format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    filemode='a',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

pep_dict = SeqIO.to_dict(SeqIO.parse(input_pep,"fasta"))
vep_df = pd.read_csv(input_vep, sep="\t",dtype={'AMINO_POS': str})
missense_df = vep_df[vep_df['VARIANT_TYPE'] == 'NONSYNONYMOUS']
missense_df['HGVS'] = missense_df['REF_AMINO'] + missense_df['AMINO_POS'] + missense_df['ALT_AMINO']
with open(sp_name + "provean_out.txt","w",buffering=1) as provean_out:
    for gene_name, group in missense_df.groupby('TRANSCRIPT_ID'):
        # if gene_name not in finished_set:
        # write fasta file for propean input
        with open(sp_name + "_provean_input.fasta.temp", "w") as pep_fasta:
            pep_fasta.write(f'>{gene_name}\n{pep_dict[gene_name].seq}\n'.replace('*', ''))
        # prepare variant file 
        variants = group['HGVS'].tolist()
        variants.sort(key=lambda z: int(z[1:-1]))
        with open(sp_name + "_provean_input.var.temp","w") as var_file:
            for variant in variants:
                var_file.write(variant + "\n")
        #run provean
        if os.path.exists("Provean_supporting_set/{gene_name}.sss"):
            provean_cmd =f'{provean_path} -q {sp_name}_provean_input.fasta.temp -v {sp_name}_provean_input.var.temp --supporting_set Provean_supporting_set/{gene_name}.sss --num_threads 2'
        else:
            provean_cmd =f'{provean_path} -q {sp_name}_provean_input.fasta.temp -v {sp_name}_provean_input.var.temp --save_supporting_set Provean_supporting_set/{gene_name}.sss --num_threads 2'
        provean_output = os.popen(provean_cmd).readlines()
        try:
            header_index = provean_output.index('# VARIATION\tSCORE\n')
            provean_scores = [x.strip().split() for x in provean_output[header_index+1:]]
            for var_info in provean_scores:
                var_name = var_info[0]
                provean_score = float(var_info[1])
                provean_out.write(f'{gene_name}\t{var_name}\t{provean_score}\n')
            logger.info(f'{gene_name} finished')
        except ValueError:
            logger.info(f'{gene_name} failed')
            logger.error(traceback.format_exc())
