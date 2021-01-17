################Cirrhosis BBJ################
rm(list = ls())
gc()
library(bigreadr)
library(dplyr)

data_path="/home/integ_project/rolypoly/"
GWAS_path<-(paste0(data_path,"GWAS_data/Cirrhosis_EAS/"))


#GWAS data for Rolypoly
GWAS_data<-fread2(paste0(GWAS_path,"Cirrhosis.auto.rsq07.mac10.txt"))
GWAS_data$MAF<-GWAS_data$MAC/(GWAS_data$N*2)
GWAS_data<-GWAS_data[,c(1,2,3,9,10,21)]
colnames(GWAS_data)<-c("chrom","pos","rsid","beta","se","maf")
GWAS_data<-subset(GWAS_data,GWAS_data$maf>0.01)
fwrite2(GWAS_data,file = paste0(GWAS_path,"GWAS_data.txt"))

#output with rsID for LDSC processing
rs_bim <- fread2("/home/comparisonProject/EAS/merge.bim")
dd <- merge(GWAS_data, rs_bim, by.x = c("CHR", "POS"), by.y = c("V1", "V4"))

dd$rsID<-dd$V2
rs_GWAS<-dd[,c("CHR","POS","rsID","Allele1","Allele2","N","BETA","SE","p.value","MAF")]
write.table(rs_GWAS,file = "/home/integ_project/rolypoly/GWAS_data/Cirrhosis_EAS/rs_GWAS.txt",sep="\t",row.names = F,quote = F)

################HCC BBJ################
rm(list = ls())
gc()
library(bigreadr)
library(dplyr)

data_path="/home/integ_project/rolypoly/"
GWAS_path<-(paste0(data_path,"GWAS_data/HCC/"))


#GWAS data for Rolypoly
GWAS_data<-fread2(paste0(GWAS_path,"HepCa.auto.rsq07.mac10.txt"))
GWAS_data$MAF<-GWAS_data$MAC/(GWAS_data$N*2)
GWAS_data<-GWAS_data[,c(1,2,3,9,10,21)]
colnames(GWAS_data)<-c("chrom","pos","rsid","beta","se","maf")
GWAS_data<-subset(GWAS_data,GWAS_data$maf>0.01)
fwrite2(GWAS_data,file = paste0(GWAS_path,"GWAS_data.txt"))

#output with rsID for LDSC processing
rs_bim <- fread2("/home/comparisonProject/EAS/merge.bim")
dd <- merge(GWAS_data, rs_bim, by.x = c("CHR", "POS"), by.y = c("V1", "V4"))

dd$rsID<-dd$V2
rs_GWAS<-dd[,c("CHR","POS","rsID","Allele1","Allele2","N","BETA","SE","p.value","MAF")]
write.table(rs_GWAS,file = "/home/integ_project/rolypoly/GWAS_data/HCC/rs_GWAS.txt",sep="\t",row.names = F,quote = F)

################Cirrhosis Gene Atlas(GA################
library(bigreadr)
library(dplyr)
library(data.table)

data_path="/home/integ_project/rolypoly/"
GWAS_path<-(paste0(data_path,"GWAS_data/Cirrhosis_EUR/"))

GA_Cirrhosis<-data.frame()
for (chr in 1:22) {
  clinical_info<-fread2(paste0(GWAS_path,"clinical_c_K74/imputed.allWhites.clinical_c_K74.chr",chr,".csv"))
  snp_info<-fread2(paste0(GWAS_path,"snp_info/snps.imputed.chr",chr,".csv"))
  snp_data<-merge(snp_info,clinical_info,by="SNP")
  snp_data$CHR<-chr
  GA_Cirrhosis<-rbind(GA_Cirrhosis,snp_data)
}

GA_Cirrhosis<-GA_Cirrhosis[,c(13,1,2,3,4,5,6,7,10,11,12)]
colnames(GA_Cirrhosis)[9:11]<-c("BETA","SE","P")
fwrite2(GA_Cirrhosis,file = paste0(GWAS_path,"GA_Cirrhosis.txt"),sep="\t")

#output for LDSC processing
rs_GWAS<-subset(GA_Cirrhosis,GA_Cirrhosis[,6]>=0.01&GA_Cirrhosis[,7]>=1E-6)
fwrite2(rs_GWAS,file = paste0(GWAS_path,"rs_GWAS.txt"),sep="\t")

#GWAS data for Rolypoly
GWAS_data<-rs_GWAS[,c(1,3,2,9,10,6)]
colnames(GWAS_data)<-c("chrom","pos","rsid","beta","se","maf")
fwrite2(GWAS_data,file = paste0(GWAS_path,"GWAS_data.txt"))
