###data preparation
#GWAS data load and filter
HepCa_GWAS<-read.table("GWAS data/HCC/HepCa.auto.rsq07.mac10.txt",sep = " ",header = T,
                       colClasses=c("numeric","numeric","character",rep("NULL",3),"numeric","NULL" ,"numeric","numeric","NULL","numeric",rep("NULL",8)))%>%funp
gwas_filter<-HepCa_GWAS[HepCa_GWAS$p.value<0.05,]
GWAS_data<-gwas_filter2[,c(1,2,3,5,6,4)]
colnames(GWAS_data)<-c("chrom","pos","rsid","beta","se","maf")
GWAS_data$maf[GWAS_data$maf>0.5]<-1-GWAS_data$maf[GWAS_data$maf>0.5]
save(GWAS_data,file = "GWAS data/HCC/GWAS_data.RData")

#single-cell data load and filter
require(Matrix)
require(dplyr)
mat_dir<-"single cell data/HCC/GSE140228/"

barcode.path <- paste0(mat_dir, 'barcodes.tsv')
gene.path <- paste0(mat_dir, 'genes.tsv')
cellinfo.path <- paste0(mat_dir, 'cellinfo.tsv')
matrix.path <- paste0(mat_dir, 'matrix.mtx')
gene.feature = read.delim(gene.path,
                          header = T,
                          stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
mat <- readMM(file = matrix.path)
colnames(mat) = barcode.names$V1
rownames(mat) = gene.feature$SYMBOL

#block annotation 
library("biomaRt")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list<-gene.feature$SYMBOL
gene_bp<-getBM(attributes = c("chromosome_name", "start_position","end_position", "external_gene_name","transcription_start_site"), 
                    filters = "external_gene_name", values = gene_list, mart = ensembl)
c<-c(1:22)
c<-as.character(c)
gene_bp_filter<-gene_bp[which(gene_bp$chromosome_name %in% c),]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$external_gene_name), ]

miRNA_idx<-grep("^MIR",gene_bp_filter$external_gene_name)
lncRNA_idx<-grep("^LINC",gene_bp_filter$external_gene_name)
uk_idx<-grep(".(\\.|\\-|\\_).",gene_bp_filter$external_gene_name)
rm_idx<-unique(c(miRNA_idx,lncRNA_idx,uk_idx))

gene_bp_filter<-gene_bp_filter[-c(rm_idx),]
block_annotation<-gene_bp_filter

block_annotation$start_position<-block_annotation$transcription_start_site-5000
block_annotation$end_position<-block_annotation$transcription_start_site+5000
block_annotation<-block_annotation[,1:4]
colnames(block_annotation)<-c("chrom","start","end","label")

save(block_annotation,file =paste0(mat_dir,"block_annotation.RData"))
load(paste0(mat_dir,"block_annotation.RData"))

#gene and cell filter
mat_f<-mat[block_annotation$label,]

cellinfo = read.delim(cellinfo.path,
                      header = T,
                      stringsAsFactors = FALSE)
type_sub<-strsplit(as.vector(cellinfo$celltype_sub),"-",fixed = F)
cell_type<-vector()
for (i in 1:length(type_sub)) {
  cell_type[i]<-type_sub[[i]][1]
}
cellinfo$cell_type<-cell_type
cellinfo$cell_type[cellinfo$celltype_sub=="Lymphoid-B-Plasma"]<-"Lymphoid-BP"
cellinfo_f<-subset(cellinfo,cellinfo$Tissue %in% c("Tumor","Normal"))
cellinfo_f<-cellinfo_f[,c("Barcode","cell_type")]
mat_f<-mat_f[,cellinfo_f$Barcode]

nzero<-vector(mode = "integer")
nzero<-apply(mat_f,1,function(x){sum(x==0)})
ncell<-nrow(mat_f)
mat_f<-mat_f[which(nzero<=0.999*ncell),]

exp_mat<-t(mat_f)%>%as.matrix

##quantile normalization:
library(preprocessCore)
matrix_expression_data<-exp_mat
normalized_exp<-normalize.quantiles(matrix_expression_data,copy=TRUE)
normalized_exp<-data.frame(normalized_exp)
colnames(normalized_exp)<-rownames(mat_f)
rownames(normalized_exp)<-colnames(mat_f)

normalized_exp$Barcode<-rownames(normalized_exp)
normalized_exp<-merge(cellinfo_f,normalized_exp,by="Barcode")
normalized_exp<-normalized_exp[,-1]
cell_norm_exp<-aggregate(normalized_exp[,2:ncol(normalized_exp)],
                         by=list(normalized_exp$cell_type),
                         FUN = mean)
rownames(cell_norm_exp)<-cell_norm_exp[,1]
cell_norm_exp<-cell_norm_exp[,-1]
cell_norm_exp<-t(cell_norm_exp)%>%as.data.frame

#data scale
cell_ns_exp<-apply(cell_norm_exp,
                     2,
                     function(x){scale(x,center = T,scale = T)})%>%abs%>%as.data.frame
rownames(cell_ns_exp)<-rownames(cell_norm_exp)
save(cell_ns_exp,file=paste0(mat_dir,"cell_ns_exp.RData"))

###Rolling rolypoly
##load data
load("GWAS data/HCC/GWAS_data.RData")
ld_path <- "EAS_LD"
load(paste0(mat_dir,"cell_ns_exp.RData"))
load(paste0(mat_dir,"block_annotation.RData"))
block_annotation_f<-subset(block_annotation,block_annotation$label %in% rownames(cell_ns_exp))

require(rolypoly)

rp <- rolypoly_roll(
  gwas_data = GWAS_data,
  block_annotation = block_annotation_f,
  block_data = cell_ns_exp,
  ld_folder = ld_path
)

rp$full_results$parameters %>% sort
rp$bootstrap_results %>% arrange(-bt_value) %>% head
plot_rolypoly_annotation_estimates(rp)
plot_rolypoly_annotation_ranking(rp)

