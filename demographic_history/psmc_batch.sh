#!/bin/bash

#PBS -N psmc_Dis4
#PBS -l nodes=1:ppn=1,walltime=168:00:00
module load samtools
module load bcftools

specieslist=""
for sample in ${specieslist}; do
cd /scratch/yufeng/demo_psmc/
mkdir ${sample}
cd ${sample}
samtools mpileup -C50 -uf /scratch/yufeng/demo_psmc/Dis_genome.fa /scratch/yufeng/demo_psmc/filteredbam/${sample}filtered_dupRM.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${sample}_genome.fq.gz
/scratch/yufeng/softwares/psmc/utils/fq2psmcfa -q20 ${sample}_genome.fq.gz > ${sample}_genome.psmcfa
/scratch/yufeng/softwares/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sample}_genome.psmc ${sample}_genome.psmcfa
/scratch/yufeng/softwares/psmc/utils/psmc_plot.pl ${sample}_genome ${sample}_genome.psmc 

/scratch/yufeng/softwares/psmc/utils/splitfa ${sample}_genome.psmcfa > ${sample}_split.psmcfa
seq 100 | xargs -i echo /scratch/yufeng/softwares/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ${sample}_round-{}.psmc ${sample}_split.psmcfa | sh
cat ${sample}_genome.psmc ${sample}_round-*.psmc > ${sample}_combined.psmc
/scratch/yufeng/softwares/psmc/utils/psmc_plot.pl ${sample}_combined ${sample}_combined.psmc
done
