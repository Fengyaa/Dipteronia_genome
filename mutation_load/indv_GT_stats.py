#!/usr/bin/env python3

#indv_GT_stats.py
import sys
import vcf

input_vcf = open(sys.argv[1],"r")
vcf_reader = vcf.Reader(open(sys.argv[1],"r"))

try:
    Type = sys.argv[2] #Radical Moderately_Radical Moderately_Conservative Conservative
except IndexError:
    Type = ""

sample_list = vcf_reader.samples
# sample_list = [x.strip() for x in open(sys.argv[3],"r")]
indexdict={}
GTdict1={}
GTdict2={}
for sample in sample_list:
    GTdict1[sample] = []
    GTdict2[sample] = []


for line in input_vcf:
    if line[0:2] == "##":
        continue
    elif line[0:2]=="#C":
        a = line.strip().split("\t")
        for sample in sample_list:
            # print(a.index(sample))
            indexdict[sample]=a.index(sample)
    else:
        # print(line)
        b = line.strip().split("\t")
        affect = b[7].split("|")[-1]
        EFF = b[7].split(";")[-1]
        # print(EFF)
        # if affect == Type or Type in EFF:
        for sample in sample_list:
            GT = b[indexdict[sample]]
            GTdict1[sample].append(GT)
            # if Type in line:
            #     # print(indexdict[sample])
            #     GTdict2[sample].append(GT)


print("indv\thom0\thom1\thet\tmissing\ttotal\tcount_1")
# print("indv\thet_ratio")
for sample in sample_list:
    hom0 = GTdict1[sample].count("0/0")
    hom1 = GTdict1[sample].count("1/1")
    het = GTdict1[sample].count("0/1") + GTdict2[sample].count("1/0")
    missing = GTdict1[sample].count("\./\.")
    total= hom0+hom1+het+missing
    count_1=het +hom1*2
    # ahom1 = GTdict1[sample].count("1/1")
    # ahet = GTdict1[sample].count("0/1") + GTdict1[sample].count("1/0")
    # acount_1=ahet/2 + ahom1
    print(sample + "\t" +str(hom0) + "\t" + str(hom1) + "\t" + str(het) + "\t" + str(missing) + "\t" + str(total) + "\t" + str(count_1))
    # print(sample + "\t" +str(tcount_1) + "\t" + str(acount_1))
    