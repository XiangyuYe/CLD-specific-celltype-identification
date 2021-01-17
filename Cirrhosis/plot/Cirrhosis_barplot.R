########
library(ggplot2)
library(patchwork)
library(dplyr)
rolypoly_path="/home/integ_project/rolypoly/single_cell_data/"
ldsc_path="/home/integ_project/LDSC_test/"

###EAS
pop<-"EAS"
rolypoly_outpath<-paste0(rolypoly_path,"GSE136103_cirrhosis/cirrhosis/output/rolypoly/rolypoly_",pop,".txt")
ldsc_outpath<-paste0(ldsc_path,"GSE136103_cirrhosis/outcome/",pop,"/Cirrhosis.cell_type_results.txt")

rolypoly_output<-read.table(rolypoly_outpath,sep = "\t",header = T)[,c("annotation","bp_value")]
rolypoly_output$Method<-"Rolypoly"
ldsc_output<-read.table(ldsc_outpath,sep = "\t",header = T)[,c("Name","Coefficient_P_value")]
colnames(ldsc_output)<-c("annotation","bp_value")
ldsc_output$Method<-"LDSC-cts"
combin_output<-rbind(rolypoly_output,ldsc_output)
combin_output<-subset(combin_output,combin_output$annotation!="intercept")
combin_output[which(combin_output$annotation=="CD8.T"),]$annotation<-"CD8+T"
combin_output[which(combin_output$annotation=="CD4.T"),]$annotation<-"CD4+T"

combin_output$sig<-ifelse(combin_output$bp_value<0.05,1,0)%>%factor(ordered = F)
nr<-nrow(combin_output)/2
combin_output$annotation<-factor(combin_output$annotation,ordered = T,
                                 levels = combin_output$annotation[nr:1])

###plot p1
p<-ggplot(data=combin_output,aes(x=annotation,y=-log10(bp_value)))
p<-p + geom_bar(aes(fill=Method),stat='identity',position = position_dodge(0.8),width = 0.8,show.legend = F)+
  scale_fill_manual(values=c("#DC0000B2","#4DBBD5B2"))
p<-p + geom_hline(yintercept = -log10(0.025),color="gold",lty=5,lwd=1)
p<-p + labs(x= "cell types",y = "-log10(P value)")+
  scale_y_continuous(limits=c(0,3),breaks =c(seq(0,3,0.5)) ) +
  scale_x_discrete(breaks=combin_output$annotation,
                   labels=combin_output$annotation)+coord_flip()+
  ggtitle("EAS")
p1<-p+theme(panel.background=element_rect(fill = "white",color = "black"),
            legend.text=element_text(size = 15,color = "grey20"),
            legend.title=element_blank(),
            panel.grid=element_blank(),
            axis.title = element_text(size = 20,color = "black"),
            axis.text = element_text(size = 18,color = "grey20"),
            plot.title = element_text(size = 20, face = "bold",hjust = 0.5))


###EUR
pop<-"EUR"
rolypoly_outpath<-paste0(rolypoly_path,"GSE136103_cirrhosis/cirrhosis/output/rolypoly/rolypoly_",pop,".txt")
ldsc_outpath<-paste0(ldsc_path,"GSE136103_cirrhosis/outcome/",pop,"/Cirrhosis.cell_type_results.txt")

rolypoly_output<-read.table(rolypoly_outpath,sep = "\t",header = T)[,c("annotation","bp_value")]
rolypoly_output$Method<-"Rolypoly"
ldsc_output<-read.table(ldsc_outpath,sep = "\t",header = T)[,c("Name","Coefficient_P_value")]
colnames(ldsc_output)<-c("annotation","bp_value")
ldsc_output$Method<-"LDSC-cts"
combin_output<-rbind(rolypoly_output,ldsc_output)
combin_output<-subset(combin_output,combin_output$annotation!="intercept")
combin_output[which(combin_output$annotation=="CD8.T"),]$annotation<-"CD8+T"
combin_output[which(combin_output$annotation=="CD4.T"),]$annotation<-"CD4+T"

combin_output$sig<-ifelse(combin_output$bp_value<0.05,1,0)%>%factor(ordered = F)
nr<-nrow(combin_output)/2
combin_output$annotation<-factor(combin_output$annotation,ordered = T,
                                 levels = combin_output$annotation[nr:1])

###plot p2
p<-ggplot(data=combin_output,aes(x=annotation,y=-log10(bp_value)))
p<-p + geom_bar(aes(fill=Method),stat='identity',position = position_dodge(0.8),width = 0.8,show.legend = T)+
  scale_fill_manual(values=c("#DC0000B2","#4DBBD5B2"))
p<-p + geom_hline(yintercept = -log10(0.025),color="gold",lty=5,lwd=1)
p<-p + labs(x= "cell types",y = "-log10(P value)")+
  scale_y_continuous(limits=c(0,3),breaks =c(seq(0,3,0.5)) ) +
  scale_x_discrete(breaks=combin_output$annotation,
                   labels=combin_output$annotation)+coord_flip()+
  ggtitle("EUR")
p2<-p+theme(panel.background=element_rect(fill = "white",color = "black"),
            legend.text=element_text(size = 15,color = "grey20"),
            legend.title=element_blank(),
            panel.grid=element_blank(),
            axis.title.y=element_blank(),
            axis.title = element_text(size = 20,color = "black"),
            axis.text = element_text(size = 18,color = "grey20") ,
            plot.title = element_text(size = 20,face = "bold",hjust = 0.5))

###plot out
pdf(paste0("/home/integ_project/rolypoly/single_cell_data/GSE136103_cirrhosis/estimates_cir.pdf"),height = 9,width = 12)
p1+p2
dev.off()

tiff(paste0("/home/integ_project/rolypoly/single_cell_data/GSE136103_cirrhosis/estimates_cir.tiff"),
     height = 9,width = 12, units = "in", res = 300,compression="lzw")
p1+p2
dev.off()
