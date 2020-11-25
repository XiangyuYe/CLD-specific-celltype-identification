rm(list=ls())
gc()
library(rolypoly)
library(dplyr)
library(ggplot2)
library(patchwork)
data_path="/net/mulan/disk2/yasheng/test/rolypoly/"

####rolling rolypoly#####
GWAS_path<-(paste0(data_path,"GWAS_data/HCC/"))
GWAS_data<-fread2(paste0(GWAS_path,"GWAS_data.txt"))
ld_path <-(paste0(data_path,"EAS_LD"))

sc_path1<-(paste0(data_path,"single_cell_data/GSE149614_HCC/"))
sc_path2<-(paste0(data_path,"single_cell_data/GSE112271_multiple/"))
sc_path3<-(paste0(data_path,"single_cell_data/integrate_HCC/"))
sc_list<-list(sc_path1,sc_path2,sc_path3)

for (i in 1:length(sc_list)) {
  
  load(paste0(sc_list[[i]],"output/rolypoly/block_annotation_all.RData"))
  load(paste0(sc_list[[i]],"output/rolypoly/cell_exp.RData"))
  
  rolypoly_gene<-intersect(rownames(cell_exp),block_annotation$label)
  miRNA_idx<-grep("^MIR",rolypoly_gene)
  lncRNA_idx<-grep("^LINC",rolypoly_gene)
  uk_idx<-grep(".(\\.|\\-|\\_).",rolypoly_gene)
  rm_idx<-unique(c(miRNA_idx,lncRNA_idx,uk_idx))
  rolypoly_gene<-rolypoly_gene[-c(rm_idx)]
  
  cell_exp_f<-cell_exp[rolypoly_gene,]
  block_annotation_f<-subset(block_annotation,block_annotation$label %in% rolypoly_gene)
  
  ###
  rp <- rolypoly_roll(
    gwas_data = GWAS_data,
    block_annotation = block_annotation_f,
    block_data = cell_exp_f,
    ld_folder = ld_path,
    bootstrap_iters =1000
  )
  
  ###outcomes extraction
  rolypoly_output<-rp$bootstrap_results
  rolypoly_output$parameters<-rp$full_results$parameters
  rolypoly_output<-rolypoly_output[order(rolypoly_output$bt_value,decreasing=T),]
  ###rolypoly_output
  write.table(rolypoly_output,file=paste0(sc_list[[i]],"output/rolypoly/test/rolypoly.txt"),
              sep = '\t',row.names = F,quote = F)
  
  rolypoly_output<-subset(rolypoly_output,rolypoly_output$annotation!="intercept")
  rolypoly_output$sig<-ifelse(rolypoly_output$CI_lo>0,1,0)%>%factor(ordered = F)
  rolypoly_output$annotation<-factor(rolypoly_output$annotation,ordered = T,levels = rolypoly_output$annotation[order(rolypoly_output$bt_value)])
  
  rolypoly_output<-subset(rolypoly_output,rolypoly_output$annotation!="intercept")
  rolypoly_output$sig<-ifelse(rolypoly_output$CI_lo>0,1,0)%>%factor(ordered = F)
  rolypoly_output$annotation<-factor(rolypoly_output$annotation,ordered = T,levels = rolypoly_output$annotation[order(rolypoly_output$bt_value)])
  
  p<-ggplot(rolypoly_output,aes(x=annotation,y=bootstrap_estimate))
  p<-p + geom_linerange(aes(ymin = CI_lo, ymax = CI_hi),lwd=2,color="grey70")
  p<-p + geom_point(aes(color=sig),size=5,show.legend = F)+
    scale_color_manual(values=c("1" = "tomato", "0" = "grey50"))
  p<-p + geom_hline(yintercept = 0,color="black",lty=5,lwd=1)
  p<-p + labs(x= "cell types",y = "estimation")+
    scale_y_continuous(limits=c(-0.0008,0.0008),n.breaks = 9)+
    coord_flip()
  p_list[[i]]<-p+theme(panel.background=element_rect(fill = "white",color = "black",size=1),
                       panel.grid=element_blank(),
                       axis.title = element_text(size = 20,color = "black"),
                       axis.text = element_text(size = 15,color = "grey20"))
}

###picture out
pdf(paste0(sc_path1,"output/rolypoly/test/estimates.pdf"),height = 12,width = 7)
  pic<-p_list[[1]]+p_list[[2]]+p_list[[3]]+plot_annotation(tag_levels = "A")
  print(pic)
dev.off()
