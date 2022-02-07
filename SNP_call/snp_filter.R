##Dis
require(data.table)
snp = fread("Dis_all.raw.vcf.snp.stat")
keep.snp <- snp$QD>=5 & !is.na(snp$QD) & snp$MQ>=40 & !is.na(snp$MQ) & snp$FS<10 & !is.na(snp$FS) & snp$SOR < 2 & !is.na(snp$SOR) & ((snp$MQR > -3 & snp$MQR<3) | is.na(snp$MQR)) & ((snp$RPR > -2 & snp$RPR<2) | is.na(snp$RPR))
keep.snp <- keep.snp & snp$Inb>-0.4 & !is.na(snp$Inb) & snp$EH<5 & !is.na(snp$EH)
keep.snp <- keep.snp & snp$miss_rate <=0.4 & snp$mean_depth>=5 & snp$mean_depth<200 & snp$mean_GQ>=20
pos.snp <-snp[keep.snp, c('chr','pos')]
write.table(pos.snp,file='Dis_all.snps.filt.txt',quote=F,row.names=F,sep="\t")

##Did
snp = fread("Did_all.raw.vcf.snp.stat")
keep.snp <- snp$QD>=5 & !is.na(snp$QD) & snp$MQ>=40 & !is.na(snp$MQ) & snp$FS<10 & !is.na(snp$FS) & snp$SOR < 2 & !is.na(snp$SOR) & ((snp$MQR > -3 & snp$MQR<3) | is.na(snp$MQR)) & ((snp$RPR > -2 & snp$RPR<2) | is.na(snp$RPR))
keep.snp <- keep.snp & snp$Inb>-0.4 & !is.na(snp$Inb) & snp$EH<5 & !is.na(snp$EH)
keep.snp <- keep.snp & snp$miss_rate <=0.4 & snp$mean_depth>=5 & snp$mean_depth<200 & snp$mean_GQ>=20
pos.snp <-snp[keep.snp, c('chr','pos')]
write.table(pos.snp,file='Did_all.snps.filt.txt',quote=F,row.names=F,sep="\t")
