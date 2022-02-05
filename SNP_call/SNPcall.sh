#!/bin/bash

#PBS -N Did_recall
#PBS -l nodes=1:ppn=4,walltime=168:00:00,mem=200g

module load bowtie2
module load samtools
module load jre1.8.0_40
cd /scratch/yufeng/reseq/kmer/
#bowtie2-build Did.genome.fa Did.genome
cat trimlist|parallel "java -jar /scratch/yufeng/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 reads/{}_1_clean.fq.gz reads/{}_2_clean.fq.gz reads/{}.1_trimmed.fq.gz reads/{}.1_unpaired.fq.gz reads/{}.2_trimmed.fq.gz reads/{}.2_unpaired.fq.gz ILLUMINACLIP:/scratch/yufeng/softwares/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && rm reads/{}_1_clean.fq.gz reads/{}_2_clean.fq.gz reads/{}.1_unpaired.fq.gz reads/{}.2_unpaired.fq.gz && bowtie2 -x Did.genome -1 reads/{}.1_trimmed.fq.gz -2 reads/{}.2_trimmed.fq.gz > alignments/{}.sam && samtools view -Sb alignments/{}.sam| samtools sort - > alignments/{}.bam && rm alignments/{}.sam"

cd alignments
##add bam tags
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=BGC07_L2.bam O=BGC07_R.bam RGID=BGC07_L2 RGLB=BGC07 RGPL=illumina RGPU=L2 RGSM=BGC07
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=BGC13_L2.bam O=BGC13_R.bam RGID=BGC13_L2 RGLB=BGC13 RGPL=illumina RGPU=L2 RGSM=BGC13
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=CJZ06_L3.bam O=CJZ06_R.bam RGID=CJZ06_L3 RGLB=CJZ06 RGPL=illumina RGPU=L3 RGSM=CJZ06
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=CJZ14_L3.bam O=CJZ14_R.bam RGID=CJZ14_L3 RGLB=CJZ14 RGPL=illumina RGPU=L3 RGSM=CJZ14
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=CJZ16_L1.bam O=CJZ16_R.bam RGID=CJZ16_L1 RGLB=CJZ16 RGPL=illumina RGPU=L1 RGSM=CJZ16
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=CPT13_L1.bam O=CPT13_R.bam RGID=CPT13_L1 RGLB=CPT13 RGPL=illumina RGPU=L1 RGSM=CPT13
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=HTS20_L4.bam O=HTS20_R.bam RGID=HTS20_L4 RGLB=HTS20 RGPL=illumina RGPU=L4 RGSM=HTS20
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=HTS23_L3.bam O=HTS23_R.bam RGID=HTS23_L3 RGLB=HTS23 RGPL=illumina RGPU=L3 RGSM=HTS23
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ09_L1.bam O=LZ09_R.bam RGID=LZ09_L1 RGLB=LZ09 RGPL=illumina RGPU=L1 RGSM=LZ09
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ12_L3.bam O=LZ12_L3_R.bam RGID=LZ12_L3 RGLB=LZ12 RGPL=illumina RGPU=L3 RGSM=LZ12
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ12_L4.bam O=LZ12_L4_R.bam RGID=LZ12_L4 RGLB=LZ12 RGPL=illumina RGPU=L4 RGSM=LZ12
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ14_L1.bam O=LZ14_R.bam RGID=LZ14_L1 RGLB=LZ14 RGPL=illumina RGPU=L1 RGSM=LZ14
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ22_L3.bam O=LZ22_R.bam RGID=LZ22_L3 RGLB=LZ22 RGPL=illumina RGPU=L3 RGSM=LZ22
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ25_L3.bam O=LZ25_R.bam RGID=LZ25_L3 RGLB=LZ25 RGPL=illumina RGPU=L3 RGSM=LZ25
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ26_L1.bam O=LZ26_R.bam RGID=LZ26_L1 RGLB=LZ26 RGPL=illumina RGPU=L1 RGSM=LZ26
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=LZ31_L3.bam O=LZ31_R.bam RGID=LZ31_L3 RGLB=LZ31 RGPL=illumina RGPU=L3 RGSM=LZ31
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=MLP01_L2.bam O=MLP01_L2_R.bam RGID=MLP01_L2 RGLB=MLP01 RGPL=illumina RGPU=L2 RGSM=MLP01
java -jar -Xmx4g picard.jar AddOrReplaceReadGroups I=MLP01_L3.bam O=MLP01_L3_R.bam RGID=MLP01_L3 RGLB=MLP01 RGPL=illumina RGPU=L3 RGSM=MLP01
##merge bam files
java -jar -Xmx4g picard.jar MergeSamFiles I=LZ12_L3_R.bam I=LZ12_L4_R.bam O=LZ12_R.bam
java -jar -Xmx4g picard.jar MergeSamFiles I=MLP01_L2_R.bam I=MLP01_L3_R.bam O=MLP01_R.bam
cd ..
cat samplelist|parallel "samtools flagstat alignments/{}.bam >alignments/{}.bam.flagstat && samtools view -q 20 -f 0x0002 -F 0X0004 -F 0X0008 -b alignments/{}_R.bam > alignments/{}_flt.bam"
cat samplelist| parallel "java -jar /scratch/yufeng/softwares/picard.jar MarkDuplicates I=alignments/{}_flt.bam O=alignments/{}_flt_MD.bam M=alignments/{}_flt_MD.log && rm alignments/{}_R.bam"

#java -jar alignments/picard.jar CreateSequenceDictionary R=Did.genome.fa O=Did.genome.dict
#samtools faidx Did.genome.fa
cat samplelist| parallel "samtools index alignments/{}_flt_MD.bam"
cat samplelist| parallel "java -Xms4g -Xmx20g -jar /export/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R Did.genome.fa -I alignments/{}_flt_MD.bam -ERC GVCF -o ./gvcfs/{}_flt_MD.g.vcf.gz"
cd gvcfs
java -Xms4g -Xmx20g -jar /export/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T CombineGVCFs -R ../Did.genome.fa -o Did_all.g.vcf.gz -V BGC07_flt_MD.g.vcf.gz -V CPT16_flt_MD.g.vcf.gz -V LZ15_flt_MD.g.vcf.gz -V MLP03_flt_MD.g.vcf.gz -V BGC13_flt_MD.g.vcf.gz -V HTS20_flt_MD.g.vcf.gz -V LZ22_flt_MD.g.vcf.gz -V  MLP04_flt_MD.g.vcf.gz -V CJZ06_flt_MD.g.vcf.gz -V HTS23_flt_MD.g.vcf.gz -V LZ25_flt_MD.g.vcf.gz -V MLP05_flt_MD.g.vcf.gz -V CJZ08_flt_MD.g.vcf.gz -V LZ06_flt_MD.g.vcf.gz -V LZ26_flt_MD.g.vcf.gz -V XMF11_flt_MD.g.vcf.gz -V CJZ14_flt_MD.g.vcf.gz -V LZ09_flt_MD.g.vcf.gz -V LZ31_flt_MD.g.vcf.gz -V CJZ16_flt_MD.g.vcf.gz -V LZ12_flt_MD.g.vcf.gz -V MLP01_flt_MD.g.vcf.gz -V CPT13_flt_MD.g.vcf.gz -V LZ14_flt_MD.g.vcf.gz -V MLP02_flt_MD.g.vcf.gz

cd /scratch/yufeng/reseq/Dis_recall
bowtie2-build Dis.genome.fa Dis.genome
cat Dis_list1|parallel "java -jar /scratch/yufeng/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 reads/{}_1_clean.fq.gz reads/{}_2_clean.fq.gz reads/{}.1_trimmed.fq.gz reads/{}.1_unpaired.fq.gz reads/{}.2_trimmed.fq.gz reads/{}.2_unpaired.fq.gz ILLUMINACLIP:/scratch/yufeng/softwares/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && rm reads/{}_1_clean.fq.gz reads/{}_2_clean.fq.gz reads/{}.1_unpaired.fq.gz reads/{}.2_unpaired.fq.gz && bowtie2 -x Dis.genome -1 reads/{}.1_trimmed.fq.gz -2 reads/{}.2_trimmed.fq.gz > alignments/{}.sam && samtools view -Sb alignments/{}.sam| samtools sort - > alignments/{}.bam && rm alignments/{}.sam"

cd alignments
##add bam tags
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM01_L4.bam O=BTM01_R.bam RGID=BTM01_L4 RGLB=BTM01 RGPL=illumina RGPU=L4 RGSM=BTM01
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM05_L2.bam O=BTM05_R.bam RGID=BTM05_L2 RGLB=BTM05 RGPL=illumina RGPU=L2 RGSM=BTM05
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM06_L4.bam O=BTM06_R.bam RGID=BTM06_L4 RGLB=BTM06 RGPL=illumina RGPU=L4 RGSM=BTM06
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM08_L3.bam O=BTM08_L3_R.bam RGID=BTM08_L3 RGLB=BTM08 RGPL=illumina RGPU=L3 RGSM=BTM08
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM08_L4.bam O=BTM08_L4_R.bam RGID=BTM08_L4 RGLB=BTM08 RGPL=illumina RGPU=L4 RGSM=BTM08
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BTM12_L3.bam O=BTM12_R.bam RGID=BTM12_L3 RGLB=BTM12 RGPL=illumina RGPU=L3 RGSM=BTM12
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BX02_L3.bam O=BX02_L3_R.bam RGID=BX02_L3 RGLB=BX02 RGPL=illumina RGPU=L3 RGSM=BX02
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BX02_L4.bam O=BX02_L4_R.bam RGID=BX02_L4 RGLB=BX02 RGPL=illumina RGPU=L4 RGSM=BX02
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BX09_L3.bam O=BX09_R.bam RGID=BX09_L3 RGLB=BX09 RGPL=illumina RGPU=L3 RGSM=BX09
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BX11_L4.bam O=BX11_R.bam RGID=BX11_L4 RGLB=BX11 RGPL=illumina RGPU=L4 RGSM=BX11
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=BX12_L2.bam O=BX12_R.bam RGID=BX12_L2 RGLB=BX12 RGPL=illumina RGPU=L2 RGSM=BX12
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=DLL07_L2.bam O=DLL07_R.bam RGID=DLL07_L2 RGLB=DLL07 RGPL=illumina RGPU=L2 RGSM=DLL07
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=DLL12_L3.bam O=DLL12_R.bam RGID=DLL12_L3 RGLB=DLL12 RGPL=illumina RGPU=L3 RGSM=DLL12
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar AddOrReplaceReadGroups I=DLL16_L2.bam O=DLL16_R.bam RGID=DLL16_L2 RGLB=DLL16 RGPL=illumina RGPU=L2 RGSM=DLL16

##merge bam files
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar MergeSamFiles I=BTM08_L3_R.bam I=BTM08_L4_R.bam O=BTM08_R.bam
java -jar -Xmx4g /scratch/yufeng/softwares/picard.jar MergeSamFiles I=BX02_L3_R.bam I=BX02_L4_R.bam O=BX02_R.bam
cd ..
cat Dis_samplelist1|parallel "samtools flagstat alignments/{}_R.bam >alignments/{}_R.bam.flagstat && samtools view -q 20 -f 0x0002 -F 0X0004 -F 0X0008 -b alignments/{}_R.bam > alignments/{}_flt.bam"
cat Dis_samplelist1|parallel "java -jar /scratch/yufeng/softwares/picard.jar MarkDuplicates I=alignments/{}_flt.bam O=alignments/{}_flt_MD.bam M=alignments/{}_flt_MD.log && rm alignments/{}_R.bam"
java -jar /scratch/yufeng/softwares/picard.jar CreateSequenceDictionary R=Dis.genome.fa O=Dis.genome.dict
samtools faidx Dis.genome.fa
cat Dis_samplelist1| parallel "samtools index alignments/{}_flt_MD.bam"
cat Dis_samplelist1| parallel "java -Xms4g -Xmx20g -jar /export/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R Did.genome.fa -I alignments/{}_flt_MD.bam -ERC GVCF -o ./gvcfs/{}_flt_MD.g.vcf.gz"
