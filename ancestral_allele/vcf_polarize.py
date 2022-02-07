#!/usr/bin/env python3
#usage python vcf_polarize.py Dis.estsfs.recode.vcf Dis_MLP05_DYD_estsfs_nomissing_out_pvalues.txt 
import sys

in_vcf = open(sys.argv[1],"r")
out_vcf = open(sys.argv[1].replace("vcf","polarized.vcf"),"w")
estsfs_out = open(sys.argv[2],"r")
ancestra_out = open(sys.argv[1].replace("vcf","ancestral.txt"),"w")
ancestra_out.write("CHROM\tPOS\tref_allele\tancestral_allele\n")


#the third colum of the estsfs output file represents the probability of the major allele being ancestral.
# 0 Site Code P-major-ancestral P-trees[A,C,G,T]
# 1 0 0.005733 0.000000 1.000000 0.000000 0.000000 
# 2 1 0.999725 0.000000 0.000000 0.000000 1.000000 
# 3 2 0.798716 0.000000 1.000000 0.000000 0.000000 

allele_list = ["A","C","G","T"]

ancestral_pval = []
for line in estsfs_out:
    if line[0] != "0":
        a = line.strip().split(" ")
        ancestral_pval.append(a[2])
estsfs_out.close()

i = -1
for line in in_vcf:
    if line[0] == "#":
        out_vcf.write(line)
    else:
        i += 1
        x = line.strip("\n").split("\t")
        CHROM = x[0]
        POS = x[1]
        ref_allele = x[3]
        alt_allele = x[4]
        allele_dic = {"A":0, "G":0, "C":0, "T":0}
        for a in range(9,len(x)):
            ingroup_P = x[a]
            ingroup_genotype = ingroup_P[0:3].split("/")
            ingroup_0_count = ingroup_genotype.count("0")
            ingroup_1_count = ingroup_genotype.count("1")                
            allele_dic[ref_allele] += ingroup_0_count
            allele_dic[alt_allele] += ingroup_1_count
        #get major allele
        major_allele = max(allele_dic, key=lambda k: allele_dic[k])
        #get minor allele 
        allele_dic.pop(major_allele)
        minor_allele = max(allele_dic, key=lambda k: allele_dic[k])
        # print(major_allele)
        # print(minor_allele)
        # print(ancestral_pval[i])
        # break
        if float(ancestral_pval[i]) > 0.9:
            ancestal_allel = major_allele
        else:
            ancestal_allel = minor_allele
        # print(ancestal_allel)
        # print(ref_allele)
        # break
        if ancestal_allel == ref_allele:
            y = [a[0:3] for a in x[9:]]
        else:
            x[3] = ancestal_allel
            x[4] = ref_allele
            trans = str.maketrans("0/1-0/0-1/1", "1/0 1/1 0/0")
            y = [a[0:3].translate(trans) for a in x[9:]]
        out_vcf.write("\t".join(x[0:9]) + "\t" + "\t".join(y) + "\n")
        ancestra_out.write(CHROM + "\t" +POS + "\t" +ref_allele + "\t" +ancestal_allel+"\n")
            

in_vcf.close()
out_vcf.close()
