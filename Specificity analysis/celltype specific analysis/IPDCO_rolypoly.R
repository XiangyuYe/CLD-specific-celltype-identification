###data path
rm(list = ls())
gc()
library(rolypoly)
library(dplyr)
library(ggplot2)
library(patchwork)

data_path="/net/mulan/disk2/yasheng/test/rolypoly/"
sc_path<-(paste0(data_path,"single_cell_data/GSE157783_IPDCO/"))

####rolling rolypoly#####
GWAS_path<-(paste0(data_path,"GWAS_data/HCC/"))
GWAS_data<-fread2(paste0(GWAS_path,"GWAS_data.txt"))

ld_path<-(paste0(data_path,"EAS_LD"))

load(paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))
load(paste0(sc_path,"output/rolypoly/cell_exp.RData"))

rolypoly_gene<-intersect(rownames(cell_exp),block_annotation$label)
miRNA_idx<-grep("^MIR",rolypoly_gene)
lncRNA_idx<-grep("^LINC",rolypoly_gene)
uk_idx<-grep(".(\\.|\\-|\\_).",rolypoly_gene)
rm_idx<-unique(c(miRNA_idx,lncRNA_idx,uk_idx))
rolypoly_gene<-rolypoly_gene[-c(rm_idx)]

cell_exp_f<-cell_exp[rolypoly_gene,]
block_annotation_f<-subset(block_annotation,block_annotation$label %in% rolypoly_gene)

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
  write.table(rolypoly_output,file=paste0(sc_path,"output/rolypoly/rolypoly_HCC.txt"),
              sep = '\t',row.names = F,quote = F)
  
