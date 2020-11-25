################################SCTransform integrate################################
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
library(scDblFinder)
options(future.globals.maxSize=4000*1024^2)
data_path="/net/mulan/disk2/yasheng/test/rolypoly/"
sc_path<-(paste0(data_path,"single_cell_data/GSE136103_cirrhosis/cirrhosis/"))

###data input and QC
filelist<-list.files(sc_path)
n<-length(filelist)-1
datalist<-list()

for (i in 1:n) {
  #i=1
  files_path<-paste0(sc_path,filelist[i])
  pre_sce<-CreateSeuratObject(counts = Read10X(files_path),min.cells = 3,
                              project = filelist[i])%>%as.SingleCellExperiment
  pdf(paste0(sc_path,"output/1.pdf"))
  pre_sce<-scDblFinder(pre_sce,verbose = T,clust.method = "overcluster")
  dev.off()
  table(pre_sce$scDblFinder.class)
  pre_seurat<-pre_sce[,pre_sce$scDblFinder.class=="singlet"]%>%as.Seurat()
  
  pre_seurat[["percent.mt"]]<-PercentageFeatureSet(pre_seurat,pattern = "^MT-")
  high.mito <- isOutlier(pre_seurat$percent.mt, type="higher")
  pre_seurat<-pre_seurat[,!high.mito]
  datalist[[i]]<-subset(pre_seurat,nFeature_RNA>=300
                        &percent.mt<=30
  )
}

my_sclist<-datalist

for (j in 1:n) {
  my_sclist[[j]]<-SCTransform(my_sclist[[j]]
                              ,variable.features.n = 3000
                              ,do.scale = F
                              ,do.center = T
                              ,return.only.var.genes = F
                              ,verbose=T)
  }

###Data Integration
my.features<-SelectIntegrationFeatures(object.list =my_sclist
                                       ,nfeatures = 3000)

my_sclist<-PrepSCTIntegration(object.list = my_sclist
                              ,anchor.features = my.features
                              ,verbose = T)

my_anchors<-FindIntegrationAnchors(object.list=my_sclist
                                   ,normalization.method = "SCT"
                                   ,anchor.features = my.features
                                   ,dims = 1:30,verbose = T,reduction = "rpca")

my_sc<-IntegrateData(anchorset=my_anchors
                     ,normalization.method = "SCT"
                     ,dims = 1:30,verbose = T)



###PCA
DefaultAssay(my_sc)<-"integrated"
my_sc<-RunPCA(my_sc,npcs = 50)

# library(scater)
# library(scran)
# percent.var <- my_sc@ reductions$ pca@ stdev
# chosen.elbow <- PCAtools::findElbowPoint(percent.var)
# chosen.elbow

pdf(file = paste0(sc_path,"output/Elbow.pdf"))
ElbowPlot(my_sc,ndims = 50)
dev.off()

my_sc<-JackStraw(my_sc,reduction = "pca",dims = 20)
my_sc<-ScoreJackStraw(my_sc,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw.pdf"))
JackStrawPlot(my_sc,dims = 1:20)
dev.off()

###clustering
n_pca=30
resolution=0.8

my_sc<-FindNeighbors(my_sc,dims=1:n_pca)
my_sc<-FindClusters(my_sc,resolution = resolution)
my_sc<-RunUMAP(my_sc,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_g",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",group.by = "orig.ident",label = T)
dev.off()

pdf(file = paste0(sc_path,"output/UMAP/UMAP",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",label = T)
dev.off()

table(my_sc@active.ident,my_sc$ident)
prop.table(table(my_sc@active.ident,my_sc$ident),2)

################################cell define################################
my_sc<-BuildClusterTree(my_sc,dims = 1:30)

pdf(file = paste0(sc_path,"output/BuildClusterTree.pdf"))
PlotClusterTree (object = my_sc)
dev.off()

###data preparation
DefaultAssay(my_sc)<-"RNA"
my_sc<-NormalizeData(my_sc,verbose = T)
my_sc<-ScaleData(my_sc,vars.to.regress = "nCount_RNA")

###marker genes
T_features_list<-c("CD2","CD3D","CD3E","CD3G")
NK_features_list<-c("KLRF1","FGFBP2","TOX2","CD3D")
MP_features_list<-c("CD14","CD163","CD68","CSF1R")
B_features_list<-c("CD79A","MS4A1","JSRP1")
Endo_features_list<-c("ENG","CDH5","VWF")
VEC_features_list<-c("RSPO3","WNT2")
LSEC_features_list<-c("CLEC4G","CLEC4M","CD14")
LEC_features_list<-c("PDPN","PROX1")
Mesenchymal_features_list<-c("COL1A2","ACTA2","COL3A1","MYH11","RGS5")
Hep_features_list<-c("ALB","FGG","APOA1","APOC3","FABP1","TM4SF4","ANXA4")
Cyc_features_list<-c("MKI67","TOP2A","UBE2C","RRM2")
MAST_features_list<-c("TPSAB1","CPA3")

names_clusters<-c("T","NK","MP","B","Endo","VEC","LSEC","LEC","Mesenchymal","Hep","cyc","MAST")
my_feature_list<-list(T_features_list,NK_features_list,MP_features_list,B_features_list,Endo_features_list,VEC_features_list,LSEC_features_list,LEC_features_list,Mesenchymal_features_list,Hep_features_list,Cyc_features_list,MAST_features_list)

##volin plot
for (m in 1:length(names_clusters)) {
  pdf(file = paste0(sc_path,"output/markers/",names_clusters[[m]],"_post.pdf"))
  VlnPlot(my_sc,features =my_feature_list[[m]]
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

##FeaturePlot
pdf(file = paste0(sc_path,"output/markers/FeaturePlot.pdf"))
FeaturePlot(my_sc,features =c("CD3D","APOC3","CD68","ENG")
        ,pt.size = 0.001,label.size = 1)%>%print()
dev.off()

##unsure cluster
my_sc.markers<-FindMarkers(my_sc,ident.1=25,only.pos = T)
head(my_sc.markers,10)

##cell type define
my_sc.new.ids<-c("T",
                 "Endothelial","NK","T","T","MP","Endothelial","Epithelial","T","B","T",
                 "MP","MP","Endothelial","MP","MP","Plasma","Cyc","Mesenchyme","Fibroblast","LEC",
                 "Epithelial","LESC","pDC","MP","T","Epithelial")
names(my_sc.new.ids)<-levels(my_sc)
my_sc<-RenameIdents(my_sc,my_sc.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_raw.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_pre<-data.frame(Barcode=colnames(my_sc),Cluster=my_sc@active.ident)
anno_pre<-subset(anno_pre,anno_pre$Cluster %in% c("Endothelial","B","Plasma","Mesenchymal","Cyc","LEC","LSEC"))
anno_pre$Cluster<-factor(anno_pre$Cluster,ordered=F)
save(anno_pre,file = paste0(sc_path,"output/anno_pre.RData"))
save(my_sc,file = paste0(sc_path,"output/rolypoly/my_sc.RData"))
###########T cells subdivide##########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_t <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="T",1]))

###
DefaultAssay(my_sc_t)<-"integrated"
my_sc_t<-RunPCA(my_sc_t,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbowt.pdf"))
ElbowPlot(my_sc_t,ndims = 50)
dev.off()

my_sc_t<-JackStraw(my_sc_t,reduction = "pca",dims = 20)
my_sc_t<-ScoreJackStraw(my_sc_t,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_t.pdf"))
JackStrawPlot(my_sc_t,dims = 1:20)
dev.off()

###
n_pca=20
resolution=0.8
my_sc_t<-FindNeighbors(my_sc_t,dims=1:n_pca)
my_sc_t<-FindClusters(my_sc_t,resolution = resolution)
my_sc_t<-RunUMAP(my_sc_t,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPt",n_pca,resolution,".pdf"))
DimPlot(my_sc_t,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_t)<-"RNA"

CD8_features_list<-c("CD8A","CD8B","GZMA","GZMH")
CD4_features_list<-c("CD4","CCR6","CCR7","SELL","FOXP3","RORC","IL17A")
NK_features_list<-c("KLRF1","FGFBP2","TOX2","CD3D")

names_T<-c("CD8","CD4","NK","cyc","Mast")
my_Tsub_list<-list(CD8_features_list,CD4_features_list,NK_features_list,Cyc_features_list,MAST_features_list)

for (m in 1:length(names_T)) {
  pdf(file = paste0(sc_path,"output/markers/",names_T[[m]],"_t.pdf"))
  VlnPlot(my_sc_t,features =my_Tsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc.markers<-FindMarkers(my_sc_t,ident.1=14,only.pos = T)
head(my_sc.markers,10)

pdf(file = paste0(sc_path,"output/markers/featureT.pdf"))
FeaturePlot(my_sc_t, c("CD4", "CD8A", "CD8B", "KLRF1","CD3D","FOXP3"))
dev.off()

my_sc_t.new.ids<-c("CD4+T",
                   "CD4+T","CD8+T","NKT","CD4+T","CD8+T","CD4+T","NKT","CD8+T","Treg","CD8+T",
                   "CD8+T","undefined","Mast","undefined")
names(my_sc_t.new.ids)<-levels(my_sc_t)
my_sc_t<-RenameIdents(my_sc_t,my_sc_t.new.ids)

DefaultAssay(my_sc_t)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_t.pdf"))
DimPlot(my_sc_t,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_t<-data.frame(Barcode=colnames(my_sc_t),Cluster=my_sc_t@active.ident)
save(anno_t,file=paste0(sc_path,"output/anno_t.RData"))
###########MP cells subdivide###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_m <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="MP",1]))
DefaultAssay(my_sc_m)

###
DefaultAssay(my_sc_m)<-"integrated"
my_sc_m<-RunPCA(my_sc_m,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_m.pdf"))
ElbowPlot(my_sc_m,ndims = 50)
dev.off()

percent.var <- my_sc_m@ reductions$ pca@ stdev
sum(percent.var[1:10])

my_sc_m<-JackStraw(my_sc_m,reduction = "pca",dims = 20)
my_sc_m<-ScoreJackStraw(my_sc_m,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_m.pdf"))
JackStrawPlot(my_sc_m,dims = 1:20)
dev.off()

###
n_pca=20
resolution=0.8
my_sc_m<-FindNeighbors(my_sc_m,dims=1:n_pca)
my_sc_m<-FindClusters(my_sc_m,resolution = resolution)
my_sc_m<-RunUMAP(my_sc_m,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPm",n_pca,resolution,".pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_m)<-"RNA"

TM_features_list<-c("S100A12","VCAN","FCN1","S100A8")
MDM_features_list<-c("MNDA","TREM2","CD9","C1QC")
DC_features_list<-c("CD1C","FCER1A","CD1E","CLEC9A")
KC_features_list<-c("CD5L","CD163","VCAM1","MARCO","CD68")

names_MP<-c("TM","MDM","DC","KC","cyc","Mast")
my_MPsub_list<-list(TM_features_list,MDM_features_list,DC_features_list,KC_features_list,
                    Cyc_features_list,MAST_features_list)

for (m in 1:length(names_MP)) {
  pdf(file = paste0(sc_path,"output/markers/",names_MP[[m]],"_m.pdf"))
  VlnPlot(my_sc_m,features =my_MPsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc_m.new.ids<-c("Monocyte",
                   "DC","Kuffer","Monocyte","MDM","Monocyte","MDM","Monocyte","Monocyte","MDM","Kuffer",
                   "pDC","CLEC9A+DC")
names(my_sc_m.new.ids)<-levels(my_sc_m)
my_sc_m<-RenameIdents(my_sc_m,my_sc_m.new.ids)

DefaultAssay(my_sc_m)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_m.pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_m<-data.frame(Barcode=colnames(my_sc_m),Cluster=my_sc_m@active.ident)
save(anno_m,file = paste0(sc_path,"output/anno_m.RData"))
###########Endo cells subdivide###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_e <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Endothelial",1]))
DefaultAssay(my_sc_e)

###
DefaultAssay(my_sc_e)<-"integrated"
my_sc_e<-RunPCA(my_sc_e,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_e.pdf"))
ElbowPlot(my_sc_e,ndims = 50)
dev.off()

my_sc_e<-JackStraw(my_sc_e,reduction = "pca",dims = 20)
my_sc_e<-ScoreJackStraw(my_sc_e,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_e.pdf"))
JackStrawPlot(my_sc_e,dims = 1:20)
dev.off()

###
n_pca=20
resolution=0.1
my_sc_e<-FindNeighbors(my_sc_e,dims=1:n_pca)
my_sc_e<-FindClusters(my_sc_e,resolution = resolution)
my_sc_e<-RunUMAP(my_sc_e,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPe",n_pca,resolution,".pdf"))
DimPlot(my_sc_e,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_e)<-"RNA"

VEC_features_list<-c("RSPO3","WNT2")
SEC_features_list<-c("CLEC4G","CLEC4M","CD14")
LEC_features_list<-c("PDPN","PROX1")

names_e<-c("VEC","LSEC","LEC")
my_Esub_list<-list(VEC_features_list,LSEC_features_list,LEC_features_list)

for (m in 1:length(names_e)) {
  pdf(file = paste0(sc_path,"output/markers/",names_e[[m]],"_e.pdf"))
  VlnPlot(my_sc_e,features =my_Esub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}


###########Epithelial cells subdivide###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_h <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Epithelial",1]))
DefaultAssay(my_sc_h)

###
DefaultAssay(my_sc_h)<-"integrated"
my_sc_h<-RunPCA(my_sc_h,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_e.pdf"))
ElbowPlot(my_sc_h,ndims = 50)
dev.off()

my_sc_h<-JackStraw(my_sc_h,reduction = "pca",dims = 20)
my_sc_h<-ScoreJackStraw(my_sc_h,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_h.pdf"))
JackStrawPlot(my_sc_h,dims = 1:20)
dev.off()

###
n_pca=20
resolution=0.4

my_sc_h<-FindNeighbors(my_sc_h,dims=1:n_pca)
my_sc_h<-FindClusters(my_sc_h,resolution = resolution)
my_sc_h<-RunUMAP(my_sc_h,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPh",n_pca,resolution,".pdf"))
DimPlot(my_sc_h,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_h)<-"RNA"
pdf(file = paste0(sc_path,"output/markers/featureH.pdf"))
FeaturePlot(my_sc_h, c("FABP1", "APOC3", "ANXA4","TM4SF4"))
dev.off()

my_sc_h.new.ids<-c("Hepatocyte",
                   "Hepatocyte","Hepatocyte","Hepatocyte","Hepatocyte","Mast")
names(my_sc_h.new.ids)<-levels(my_sc_h)
my_sc_h<-RenameIdents(my_sc_h,my_sc_h.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAPh.pdf"))
DimPlot(my_sc_h,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_h<-data.frame(Barcode=colnames(my_sc_h),Cluster=my_sc_h@active.ident)
save(anno_h,file = paste0(sc_path,"output/anno_h.RData"))

################################seurat output################################
load(paste0(sc_path,"output/anno_m.RData"))
load(paste0(sc_path,"output/anno_t.RData"))
load(paste0(sc_path,"output/anno_h.RData"))
load(paste0(sc_path,"output/anno_pre.RData"))
load(paste0(sc_path,"output/my_sc.RData"))
anno_data<-rbind(anno_pre,anno_t,anno_m,anno_h)
###filter and order
anno_data<-subset(anno_data,anno_data$Cluster!="undefined")
anno_data$Cluster<-factor(anno_data$Cluster,levels = c("Endothelial","CD4+T","Monocyte","CD8+T","Hepatocyte","NK",
                                                         "B","Macrophage","DC","Mesenchymal","Plasma","Treg",
                                                         "Cyc","LEC","pDC","LSEC","Mast"))
my_sc<-subset(my_sc,cells=as.vector(anno_data$Barcode))
anno_data<-anno_data[colnames(my_sc),]
Idents(my_sc)<-anno_data$Cluster

pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

###marker gene present
DefaultAssay(my_sc)<-"RNA"
my_sc.markers<-FindAllMarkers(my_sc,only.pos = T)
my_sc.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_logFC)

top10<-my_sc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
write.table(top10,file = paste0(sc_path,"output/post_top10_markers_f.txt"),sep = "\t")

top5<-my_sc.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap_f.pdf"),height = 10,width = 15)
DoHeatmap(my_sc,features = top5$gene,size = 2,label = F)
dev.off()

###data output
save(anno_data,file = paste0(sc_path,"output/rolypoly/anno_data.RData"))
save(my_sc,file = paste0(sc_path,"output/rolypoly/my_sc.RData"))
################################gene annotation################################
load(paste0(sc_path,"output/my_sc.RData"))
DefaultAssay(my_sc)<-"RNA"
library("biomaRt")
library(httr)
set_config(config(ssl_verifypeer = 0L))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",GRCh = 37)
ensembl_list<-read.table("/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE136103_cirrhosis/cirrhosis/1_cd45+/genes.tsv"
                         ,header = F,stringsAsFactors = F)
gene_bp<-getBM(attributes = c("ensembl_gene_id","chromosome_name", "transcript_start","transcript_end","transcription_start_site","start_position","end_position","transcript_version"), 
               filters = "ensembl_gene_id", values = ensembl_list[,1], mart = ensembl)

c<-c(1:22)
c<-as.character(c)
gene_bp_filter<-gene_bp[which(gene_bp$chromosome_name %in% c),]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$ensembl_gene_id), ]
gene_bp_filter<-gene_bp_filter[order(gene_bp_filter$transcription_start_site),]

rownames(ensembl_list)<-ensembl_list[,1]
ensembl_list<-ensembl_list[gene_bp_filter$ensembl_gene_id,]
gene_bp_filter$external_gene_name<-ensembl_list[,2]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$external_gene_name), ]
###gene_coord
gene_coord_all<-gene_bp_filter[,c("external_gene_name","chromosome_name", "start_position","end_position")]
colnames(gene_coord_all)<-c("GENE","CHR","START","END")
rownames(gene_coord_all)<-gene_coord_all$GENE
save(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/gene_coord_all.RData")
write.table(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/gene_coord/gene_coord_all.txt",
            sep = "\t",quote = F,col.names = T,row.names = F)

###block_annotation
block_annotation<-gene_bp_filter
block_annotation$start_position<-block_annotation$transcription_start_site-10000
block_annotation$end_position<-block_annotation$transcription_start_site+10000
block_annotation<-block_annotation[, c("chromosome_name", "start_position","end_position","external_gene_name")]
colnames(block_annotation)<-c("chrom","start","end","label")
save(block_annotation,file =paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))


################################rolypoly input################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/anno_data.RData"))
###gene_retain
counts_exp<-GetAssayData(my_sc,assay="RNA",slot="counts")
gene_retain<-CreateAssayObject(counts=counts_exp,min.cells = 3)%>%rownames()

###RNA matrix
RNA_exp<-GetAssayData(my_sc,assay="RNA",slot="data")[gene_retain,]%>%t()
RNA_exp<-as.data.frame(RNA_exp)

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc","Mast","CLEC9A+DC"))
anno_data$Barcode<-as.character(anno_data$Barcode)
anno_data$Cluster<-as.character(anno_data$Cluster)
exp_data<-exp_data[anno_data$Barcode,]

exp_data$Barcode<-rownames(exp_data)
exp_data<-merge(exp_data,anno_data,by="Barcode")[,-1]
#merge and aggregate with mean
exp_data<-aggregate(exp_data[,1:(ncol(exp_data)-1)],
                    list(type=exp_data$Cluster),
                    mean)
rownames(exp_data)<-exp_data[,1]
exp_data<-exp_data[,-1]
# scale
exp_data<-scale(exp_data)
exp_data<-t(exp_data)%>%as.data.frame

cell_exp<-abs(exp_data)
save(cell_exp,file = paste0(sc_path,"output/rolypoly/cell_exp.RData"))


################################LDSC_cts input################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))
my_sc_f<-my_sc[block_annotation$label,]

load(paste0(sc_path,"output/rolypoly/anno_data.RData"))
my_sc_f<- subset(my_sc_f,
                 cells=as.vector(anno_data[!anno_data$Cluster%in%c("undefined","Cyc","CLEC9A+DC","Mast"),1]))

###gene list
DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/10
DefaultAssay(my_sc_f)<-"RNA"
all_clusters<-c("Endothelial","CD4+T","Monocyte","CD8+T","Hepatocyte","NK",
                "B","Kuffer","NKT","DC","Mesenchymal","MDM","Plasma","Treg","LEC","pDC","LSEC")

gene_set<-list()
for (c in 1:length(all_clusters)) {
  ident_x<-all_clusters[c]
  ident_y<-setdiff(all_clusters,ident_x)
  gene_set[[c]]<-FindMarkers(my_sc_f,ident.1=ident_x,
                             ident.2=ident_y,
                             logfc.threshold = 0,min.pct=0.001,
                             only.pos = T)
  gene_set_f<-gene_set[[c]]%>%top_n(n=-n_top,wt=p_val)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,
              file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/gene_list/EAS/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}

gene_coord_all<-read.table(file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/gene_coord/gene_coord_all.txt",
                           sep = "\t",header = T)
###conrol list
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/gene_list/EAS/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)

###ldcts
all_clusters<-c("Endothelial","CD4+T","Monocyte","CD8+T","Hepatocyte","NK",
                "B","Kuffer","NKT","DC","Mesenchymal","MDM","Plasma","Treg","LEC","pDC","LSEC")
path_list<-(paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/ldsc_out/EAS/",
                   1:length(all_clusters),
                   "_,/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/ldsc_out/EAS/0_"))

Cirrhosis_ldcts<-data.frame(all_clusters,path_list)
write.table(Cirrhosis_ldcts,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE136103_cirrhosis/GSE136103_cirrhosis_EAS_ldcts"),
            sep = " ",quote = F,col.names =F,row.names =F)




