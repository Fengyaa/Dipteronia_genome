import sys
import os

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


anc_seq = ''
anc_pep = ''
n_base = 0
codon = ''
OG_name = sys.argv[1].split('.')[0]

with open(sys.argv[1],'r') as in_file: ## input name format :OG0017656.codon.fasta.state
    for line in in_file:
        if line[0] != '#':
            x = line.strip().split()
            if x[0] == 'Node1':
                anc_seq += x[2]
                n_base += 1
                codon += x[2]
                if len(codon) == 3:
                    if '-' in codon:
                        anc_pep += '?'
                    else:
                        anc_pep += DNA_Codons[codon]
                    codon = ''
align = read_fasta(f'{OG_name}.codon.fasta')

if len(anc_seq) != len(list(align.values())[0]):
    print(OG_name)
 
with open(f'{OG_name}.anc.cds','w') as out_cds:
    out_cds.write(anc_seq + "\n")
with open(f'{OG_name}.anc.pep','w') as out_pep:
    out_pep.write(anc_pep + "\n")
