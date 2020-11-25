###########
library(ggplot2)
library(patchwork)
library(dplyr)
rolypoly_path="/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/"
ldsc_path="/net/mulan/disk2/yasheng/test/LDSC_test/"

data_set<-"GSE157783_IPDCO"
rolypoly_outpath<-paste0(rolypoly_path,data_set,"/output/rolypoly/rolypoly_HCC.txt")
ldsc_outpath<-paste0(ldsc_path,data_set,"/outcome/EAS/HCC.cell_type_results.txt")

##data combination
rolypoly_output<-read.table(rolypoly_outpath,sep = "\t",header = T)[,c("annotation","bp_value")]
rolypoly_output$Method<-"Rolypoly"
ldsc_output<-read.table(ldsc_outpath,sep = "\t",header = T)[,c("Name","Coefficient_P_value")]
colnames(ldsc_output)<-c("annotation","bp_value")
ldsc_output$Method<-"LDSC-cts"
combin_output<-rbind(rolypoly_output,ldsc_output)
combin_output<-subset(combin_output,combin_output$annotation!="intercept")
combin_output[which(combin_output$annotation=="CADPS2.neurons"),]$annotation<-"CADPS2+neurons"

combin_output$sig<-ifelse(combin_output$bp_value<0.05,1,0)%>%factor(ordered = F)
nr<-nrow(combin_output)/2
combin_output$annotation<-factor(combin_output$annotation,ordered = T,
                                 levels = combin_output$annotation[nr:1])

###plot
p<-ggplot(data=combin_output,aes(x=annotation,y=-log10(bp_value)))
p<-p + geom_bar(aes(fill=Method),stat='identity',position = "dodge",width = 1,show.legend = T)#+
p<-p + geom_hline(yintercept = -log10(0.05),color="gold",lty=5,lwd=1)
p<-p + labs(x= "cell types",y = "-log10(P value)")+
  scale_y_continuous(limits=c(0,2),breaks =c(seq(0,2,0.5)) ) +
  scale_x_discrete(breaks=combin_output$annotation,
                   labels=combin_output$annotation)+coord_flip()
p1<-p+theme(panel.background=element_rect(fill = "white",color = "black"),
             legend.text=element_text(size = 15,color = "grey20"),
             legend.title=element_blank(),
             panel.grid=element_blank(),
             axis.title = element_text(size = 20,color = "black"),
             axis.text = element_text(size = 15,color = "grey20") )


###plot out
pdf(paste0("/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE157783_IPDCO/estimates_brain_HCC.pdf"),height = 9,width = 9)
p1
dev.off()

tiff(paste0("/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE157783_IPDCO/estimates_brain_HCC.tiff"),
     height = 9,width = 9, units = "in", res = 300,compression="lzw")
p1
dev.off()

