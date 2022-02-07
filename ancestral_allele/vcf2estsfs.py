#!/usr/bin/env python3

##usage: python input.vcf outgroupname1 outgroupname2


import sys
import re 
ingroup_vcf = open(sys.argv[1],"r")
estsfs_out = open(sys.argv[1].replace(".vcf","_estsfs.txt"),"w")
positions_out = open(sys.argv[1].replace(".vcf","_estsfs.positions.txt"),"w")
outgroupname1 = sys.argv[2]
outgroupname2 = sys.argv[3]

for line in ingroup_vcf:
    if line.startswith("#CHROM"):
        samples_list=line.strip().split()
        outgroup_index1=samples_list.index(outgroupname1)
        outgroup_index2=samples_list.index(outgroupname2)
    if line[0] != "#":
        x = line.strip("\n").split("\t")
        CHROM = x[0]
        POS = x[1]
        ref_allele = x[3]
        alt_allele = x[4]
        allele_dic = {"A":0, "C":0, "G":0, "T":0}
        if len(ref_allele) == len(alt_allele) ==1:
            x.pop(outgroup_index)
            ingroup_GT = "".join([i.split(":")[0] for i in x[9:]])
            ref_count = ingroup_GT.count("0")
            alt_count = ingroup_GT.count("1")
            allele_dic[ref_allele]=ref_count
            allele_dic[alt_allele]=alt_count
            # for a in range(9,len(x)-1):
            #     ingroup_P = x[a]
            #     ingroup_genotype_u = ingroup_P.split(":")[0]
            #     if ingroup_genotype_u == "0/0" or x == "0/1":
            #         ingroup_genotype = ingroup_genotype_u.split("/")
            #         ingroup_0_count = ingroup_genotype.count("0")
            #         ingroup_1_count = ingroup_genotype.count("1")
            #     else:
            #         ingroup_0_count = 0
            #         ingroup_1_count = 0

                # allele_dic[ref_allele] += ingroup_0_count
                # allele_dic[alt_allele] += ingroup_1_count
            ingroup = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])

            allele_dic = {"A":0, "C":0, "G":0, "T":0}
            outgroup1 = x[outgroup_index1]
            outgroup1_genotype = outgroup1[0:3].split("/")
            outgroup1_0_count = outgroup1_genotype.count("0")
            outgroup1_1_count = outgroup1_genotype.count("1")               
            allele_dic[ref_allele] = outgroup1_0_count
            allele_dic[alt_allele] = outgroup1_1_count
            #exclude het sites in outgroup
            if outgroup1_0_count ==1 or outgroup1_1_count == 1:
                allele_dic[ref_allele] = 0
                allele_dic[alt_allele] = 0
            else:
                allele_dic[ref_allele] = int(outgroup1_0_count/2)
                allele_dic[alt_allele] = int(outgroup1_1_count/2)

            outgroup1 = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])
            
            allele_dic = {"A":0, "C":0, "G":0, "T":0}
            outgroup2 = x[outgroup_index2]
            outgroup2_genotype = outgroup2[0:3].split("/")
            outgroup2_0_count = outgroup2_genotype.count("0")
            outgroup2_1_count = outgroup2_genotype.count("1")                
            allele_dic[ref_allele] = outgroup2_0_count
            allele_dic[alt_allele] = outgroup2_1_count
            if outgroup2_0_count == 1 or outgroup2_1_count == 1:
                allele_dic[ref_allele] = 0
                allele_dic[alt_allele] = 0
            else:
                allele_dic[ref_allele] = int(outgroup2_0_count/2)
                allele_dic[alt_allele] = int(outgroup2_1_count/2)
            outgroup2 = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])
            
            if outgroup1 == outgroup2 == "0,0,0,0":
                continue
            else:
                estsfs_out.write(ingroup + "\t" + outgroup1 + " " + outgroup2 + "\n")
                # estsfs_out.write(ingroup + "\t" + outgroup1 + "\n")
                positions_out.write(CHROM + "\t" + POS + "\n")



ingroup_vcf.close()
estsfs_out.close()
positions_out.close()
