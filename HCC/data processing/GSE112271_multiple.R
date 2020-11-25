#################################SCTransform And Integration#############################################################
#BiocManager::install("scater")
rm(list=ls())
gc()
library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(scDblFinder)

options(future.globals.maxSize=4000*1024^2)

data_path="/net/mulan/disk2/yasheng/test/rolypoly/"
sc_path<-(paste0(data_path,"single_cell_data/GSE112271_multiple/"))

###data input
filelist<-list.files(sc_path)
n<-length(filelist)-1
datalist<-list()

for (i in 1:n) {
  #i=1
  files_path<-paste0(sc_path,filelist[i])
  pre_sce<-CreateSeuratObject(counts = Read10X(files_path),min.cells = 3,
                              project = filelist[i])%>%as.SingleCellExperiment
  pdf(paste0(sc_path,"output/",1,".pdf"))
  pre_sce<-scDblFinder(pre_sce,verbose = T,clust.method = "overcluster")
  dev.off()
  table(pre_sce$scDblFinder.class)
  pre_seurat<-pre_sce[,pre_sce$scDblFinder.class=="singlet"]%>%as.Seurat()
  
  #Quality control
  pre_seurat[["percent.mt"]]<-PercentageFeatureSet(pre_seurat,pattern = "^MT-")
  high.mito <- isOutlier(pre_seurat$percent.mt, type="higher")
  pre_seurat<-pre_seurat[,!high.mito]

  pre_seurat<-subset(pre_seurat,nFeature_RNA>=200
                     &percent.mt<=30
  )
  pre_seurat<-SCTransform(pre_seurat
                          ,variable.features.n = 3000
                          ,vars.to.regress = c("nCount_RNA","percent.mt")
                          ,do.scale = F
                          ,do.center = T
                          ,return.only.var.genes = F
                          ,verbose=F)
  datalist[[i]]<-pre_seurat
}

###Data Integration
my.features<-SelectIntegrationFeatures(object.list =datalist
                                       ,nfeatures = 3000)

datalist<-PrepSCTIntegration(object.list = datalist
                             ,anchor.features = my.features
                             ,verbose = T)

my_anchors<-FindIntegrationAnchors(object.list=datalist
                                   ,normalization.method = "SCT"
                                   ,anchor.features = my.features
                                   ,dims = 1:30,verbose = T)

my_sc<-IntegrateData(anchorset=my_anchors
                     ,normalization.method = "SCT"
                     ,dims = 1:30,verbose = T)


###Run PCA
my_sc<-RunPCA(my_sc,npcs = 50)
pdf(file = paste0(sc_path,"output/DimPlot.pdf"))
DimPlot(my_sc,reduction = "pca",split.by  = "orig.ident",ncol = 2)
dev.off()

tiff(file = paste0(sc_path,"output/DimHeatmap.tiff"),height = 5000,3000,compression = "lzw")
DimHeatmap(my_sc,dims = c(1:9),reduction = "pca",balanced = T)
dev.off()

pdf(file = paste0(sc_path,"output/Elbow.pdf"))
ElbowPlot(my_sc,ndims = 50)
dev.off()


###clustering
n_pca=40
resolution=0.8
DefaultAssay(my_sc)<-"integrated"
###
my_sc<-FindNeighbors(my_sc,dims=1:n_pca)
my_sc<-FindClusters(my_sc,resolution = resolution,dims=1:n_pca)
my_sc<-RunUMAP(my_sc,dims = 1:n_pca)


pdf(file = paste0(sc_path,"output/UMAP/UMAP_g",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",group.by = "orig.ident",label = T)
dev.off()

pdf(file = paste0(sc_path,"output/UMAP/UMAP",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",label = T)
dev.off()



###################cell define######################
# DefaultAssay(my_sc)<-"integrated"
DefaultAssay(my_sc)<-"RNA"
my_sc<-NormalizeData(my_sc,verbose = T)
my_sc<-ScaleData(my_sc,vars.to.regress = c("nCount_RNA","percent.mt"))
#marker genes
T_features_list<-c("CD2","CD3D","CD3E","CD3G")
MP_features_list<-c("CD14","CD163","CD68","CSF1R")
B_features_list<-c("CD79A","MS4A1","JSRP1")
Endo_features_list<-c("ENG","CDH5","VWF")
Mesenchymal_features_list<-c("COL1A2","ACTA2","COL3A1","MYH11","RGS5")
Epi_features_list<-c("ALB","FGG","APOA1","APOC3","FABP1","TM4SF4","ANXA4")
Cyc_features_list<-c("MKI67","TOP2A","UBE2C","RRM2")
MAST_features_list<-c("TPSAB1","CPA3")
NK_features_list<-c("KLRF1","FGFBP2","GNLY","NKG7")

names_clusters<-c("T","MP","B","Endo","Mesenchymal","Epi","cyc","MAST","NK")
my_feature_list<-list(T_features_list,MP_features_list,B_features_list,Endo_features_list,Mesenchymal_features_list,Epi_features_list,Cyc_features_list,MAST_features_list,NK_features_list)

for (m in 1:length(names_clusters)) {
  pdf(file = paste0(sc_path,"output/markers/",names_clusters[[m]],"_post.pdf"))
  VlnPlot(my_sc,features =my_feature_list[[m]]
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/FeaturePlot.pdf"))
FeaturePlot(my_sc,features =c("CD3D","APOC3","CD68","ENG")
            ,pt.size = 0.001,label.size = 1)%>%print()
dev.off()

my_sc<-BuildClusterTree(my_sc,dims = 1:30)
pdf(file = paste0(sc_path,"output/BuildClusterTree.pdf"))
PlotClusterTree (object = my_sc)
dev.off()

my_sc.markers<-FindMarkers(my_sc,ident.1="T")
head(my_sc.markers,10)

my_sc.new.ids<-c("Epithelial",
                 "MP","Epithelial","T","Mesenchymal","Endothelial","Epithelial","Epithelial","Epithelial","Epithelial","T",
                 "Epithelial","MP","Epithelial","MP","Endothelial","Endothelial","Cyc","Fibroblast","Cyc","Mesenchymal",
                 "Cyc","B","Endothelial","MP","Endothelial")
names(my_sc.new.ids)<-levels(my_sc)
my_sc<-RenameIdents(my_sc,my_sc.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_raw.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_pre<-data.frame(Barcode=colnames(my_sc),Cluster=my_sc@active.ident)
anno_pre<-subset(anno_pre,anno_pre$Cluster %in% c("B","Cyc"))
anno_pre$Cluster<-factor(anno_pre$Cluster,ordered=F)
save(anno_pre,file =paste0(sc_path,"output/anno_pre.RData"))

save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))

###########T###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_t <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="T",1]))

###PCA
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

###Clustering
n_pca=30
resolution=0.8
DefaultAssay(my_sc_t)<-"integrated"
###
my_sc_t<-FindNeighbors(my_sc_t,dims=1:n_pca)
my_sc_t<-FindClusters(my_sc_t,resolution = resolution)
my_sc_t<-RunUMAP(my_sc_t,dims = 1:n_pca)
###
pdf(file = paste0(sc_path,"output/UMAP/UMAPt",n_pca,resolution,".pdf"))
DimPlot(my_sc_t,reduction = "umap",label = T)
dev.off()
###
DefaultAssay(my_sc_t)<-"RNA"
CD8_features_list<-c("CD8A", "CD8B","GZMA","GZMH")
CD4_features_list<-c("CD4","CCR6","CCR7","SELL","FOXP3","RORC","IL17A")
NK_features_list<-c("KLRF1","FGFBP2","TOX2","CD3D")
Cyc_features_list<-c("MKI67","TOP2A","UBE2C","RRM2")

names_T<-c("CD8","CD4","NK","cyc")
my_Tsub_list<-list(CD8_features_list,CD4_features_list,NK_features_list,Cyc_features_list)

for (m in 1:length(names_T)) {
  pdf(file = paste0(sc_path,"output/markers/",names_T[[m]],"_t.pdf"))
  VlnPlot(my_sc_t,features =my_Tsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/featureT.pdf"))
FeaturePlot(my_sc_t, c("CD4", "CD8A", "KLRF1","CD3D"))
dev.off()

my_sc_t.markers<-FindMarkers(my_sc_t,ident.1=3,only.pos = T)
head(my_sc_t.markers,10)

###########MP######
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_m <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="MP",1]))

###PCA
DefaultAssay(my_sc_m)<-"integrated"
my_sc_m<-RunPCA(my_sc_m,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_m.pdf"))
ElbowPlot(my_sc_m,ndims = 50)
dev.off()

my_sc_m<-JackStraw(my_sc_m,reduction = "pca",dims = 20)
my_sc_m<-ScoreJackStraw(my_sc_m,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_m.pdf"))
JackStrawPlot(my_sc_m,dims = 1:20)
dev.off()

###Clustering
n_pca=20
resolution=0.8
DefaultAssay(my_sc_m)<-"integrated"
###
my_sc_m<-FindNeighbors(my_sc_m,dims=1:n_pca)
my_sc_m<-FindClusters(my_sc_m,resolution = resolution)
my_sc_m<-RunUMAP(my_sc_m,dims = 1:n_pca)
###
pdf(file = paste0(sc_path,"output/UMAP/UMAPm",n_pca,resolution,".pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T)
dev.off()

#mp cells
DefaultAssay(my_sc_m)<-"RNA"

TM_features_list<-c("S100A12","VCAN","FCN1","S100A8")
MDM_features_list<-c("MNDA","TREM2","CD9","C1QC")
DC_features_list<-c("CD1C","FCER1A","CD1E","CLEC9A")
KC_features_list<-c("CD5L","CD163","VCAM1","MARCO","CD68")
Cyc_features_list<-c("MKI67","TOP2A","UBE2C","RRM2")

names_MP<-c("TM","MDM","DC","KC","cyc")
my_MPsub_list<-list(TM_features_list,MDM_features_list,DC_features_list,KC_features_list,
                    Cyc_features_list)

###
for (m in 1:length(names_MP)) {
  pdf(file = paste0(sc_path,"output/markers/",names_MP[[m]],"_m.pdf"))
  VlnPlot(my_sc_m,features =my_MPsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/featureM.pdf"))
FeaturePlot(my_sc_m, c("MARCO", "TREM2", "MNDA","CD5L"))
dev.off()

my_sc.markers<-FindMarkers(my_sc_m,ident.1=2,only.pos = T)
head(my_sc.markers,10)

my_sc_m.new.ids<-c("MDM",
                   "MDM","MDM","MDM","MDM","DC","MDM","MDM","undefined","Kuffer","MDM",
                   "Mesenchyme")
names(my_sc_m.new.ids)<-levels(my_sc_m)
my_sc_m<-RenameIdents(my_sc_m,my_sc_m.new.ids)

my_sc_m.markers<-FindAllMarkers(my_sc_m,only.pos = T)
my_sc_m.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
top10<-my_sc_m.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap_M.pdf"))
DoHeatmap(my_sc_m,features = top10_m$gene,size = 2,label=F)
dev.off()

DefaultAssay(my_sc_m)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_m.pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_m<-data.frame(Barcode=colnames(my_sc_m),Cluster=my_sc_m@active.ident)
save(anno_m,file = paste0(sc_path,"output/anno_m.RData"))

###########Endo#####
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_e <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Endothelial",1]))

###PCA
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

###Clustering
n_pca=20
resolution=0.6

my_sc_e<-FindNeighbors(my_sc_e,dims=1:n_pca)
my_sc_e<-FindClusters(my_sc_e,resolution = resolution)
my_sc_e<-RunUMAP(my_sc_e,dims = 1:n_pca)
###
pdf(file = paste0(sc_path,"output/UMAP/UMAPe",n_pca,resolution,".pdf"))
DimPlot(my_sc_e,reduction = "umap",label = T)
dev.off()

DefaultAssay(my_sc_e)<-"RNA"

VEC_features_list<-c("RSPO3","WNT2")
LSEC_features_list<-c("CLEC4G","CLEC4M","CD14")
LEC_features_list<-c("PDPN","PROX1")

names_e<-c("VEC","LSEC","LEC")
my_Esub_list<-list(VEC_features_list,LSEC_features_list,LEC_features_list)

for (m in 1:length(names_e)) {
  pdf(file = paste0(sc_path,"output/markers/",names_e[[m]],"_e.pdf"))
  VlnPlot(my_sc_e,features =my_Esub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc_e.new.ids<-c("Endothelial",
                   "Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","LEC","LSEC")
names(my_sc_e.new.ids)<-levels(my_sc_e)
my_sc_e<-RenameIdents(my_sc_e,my_sc_e.new.ids)


my_sc_e.markers<-FindAllMarkers(my_sc_e,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
top10<-my_sc_e.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap_e.pdf"))
DoHeatmap(my_sc_e,features = top10$gene,size = 2,label = F)
dev.off()

DefaultAssay(my_sc_e)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_e.pdf"))
DimPlot(my_sc_e,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_e<-data.frame(Barcode=colnames(my_sc_e),Cluster=my_sc_e@active.ident)
save(anno_e,file=paste0(sc_path,"output/anno_e.RData"))

###########H###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_h <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Epithelial",1]))

###PCA
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

###Clustering
n_pca=20
resolution=0.4
DefaultAssay(my_sc_h)<-"integrated"
###
my_sc_h<-FindNeighbors(my_sc_h,dims=1:n_pca)
my_sc_h<-FindClusters(my_sc_h,resolution = resolution)
my_sc_h<-RunUMAP(my_sc_h,dims = 1:n_pca)
###
pdf(file = paste0(sc_path,"output/UMAP/UMAPh",n_pca,resolution,".pdf"))
DimPlot(my_sc_h,reduction = "umap",label = T)
dev.off()

DefaultAssay(my_sc_h)<-"RNA"

Hep_features_list<-c("ALB","FGG","APOA1","APOC3")
Chol_features_list<-c("FABP1","TM4SF4","ANXA4")

names_h<-c("Hep","Chol")
my_Hsub_list<-list(Hep_features_list,Chol_features_list)

for (m in 1:length(names_h)) {
  pdf(file = paste0(sc_path,"output/markers/",names_h[[m]],"_h.pdf"))
  VlnPlot(my_sc_h,features =my_Hsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc.markers<-FindMarkers(my_sc_h,ident.1=4,only.pos=T)
head(my_sc.markers,10)

pdf(file = paste0(sc_path,"output/markers/featureH.pdf"))
FeaturePlot(my_sc_h, c("FABP1", "APOC3", "ANXA4","TM4SF4"))
dev.off()

my_sc_h.new.ids<-c("Hepatocyte",
                   "Hepatocyte","Hepatocyte","Hepatocyte","Cholangiocyte","Hepatocyte","Hepatocyte","Hepatocyte")
names(my_sc_h.new.ids)<-levels(my_sc_h)
my_sc_h<-RenameIdents(my_sc_h,my_sc_h.new.ids)

my_sc_h.markers<-FindAllMarkers(my_sc_h,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
top10<-my_sc_h.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap_H.pdf"))
DoHeatmap(my_sc_h,features = top10$gene,size = 2,label=F)
dev.off()

pdf(file = paste0(sc_path,"output/UMAP/UMAPh.pdf"))
DimPlot(my_sc_h,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_h<-data.frame(Barcode=colnames(my_sc_h),Cluster=my_sc_h@active.ident)
save(anno_h,file = paste0(sc_path,"output/anno_h.RData"))

###########Me###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_Me <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Mesenchymal",1]))

###PCA
DefaultAssay(my_sc_Me)<-"integrated"
my_sc_Me<-RunPCA(my_sc_Me,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_Me.pdf"))
ElbowPlot(my_sc_Me,ndims = 50)
dev.off()

my_sc_Me<-JackStraw(my_sc_Me,reduction = "pca",dims = 20)
my_sc_Me<-ScoreJackStraw(my_sc_Me,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw_Me.pdf"))
JackStrawPlot(my_sc_Me,dims = 1:20)
dev.off()

###clustering
n_pca=10
resolution=0.4
###
my_sc_Me<-FindNeighbors(my_sc_Me,dims=1:n_pca)
my_sc_Me<-FindClusters(my_sc_Me,resolution = resolution)
my_sc_Me<-RunUMAP(my_sc_Me,dims = 1:n_pca)
###
pdf(file = paste0(sc_path,"output/UMAP/UMAPMe",n_pca,resolution,".pdf"))
DimPlot(my_sc_Me,reduction = "umap",label = T)
dev.off()

DefaultAssay(my_sc_Me)<-"RNA"
Mesenchymal_features_list<-c("COL1A2","ACTA2","COL3A1","MYH11","RGS5")
VEC_features_list<-c("RSPO3","WNT2")
LSEC_features_list<-c("CLEC4G","CLEC4M","CD14")
LEC_features_list<-c("PDPN","PROX1")

names_me<-c("Mesenchymal","VEC","LSEC","LEC","")
my_mesub_list<-list(Mesenchymal_features_list,VEC_features_list,LSEC_features_list,LEC_features_list)

for (m in 1:length(names_me)) {
  pdf(file = paste0(sc_path,"output/markers/",names_me[[m]],"_Me.pdf"))
  VlnPlot(my_sc_Me,features =my_mesub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc.markers<-FindMarkers(my_sc_Me,ident.1=5,only.pos=T)
head(my_sc.markers,10)

pdf(file = paste0(sc_path,"output/markers/featureMe.pdf"))
FeaturePlot(my_sc_Me, c("COL1A1","CDH5","VWF","RAMP2"))
dev.off()

my_sc_Me.new.ids<-c("Mesenchymal",
                   "Mesenchymal","Mesenchymal","Fibroblast","Mesenchymal","Endothelial")
names(my_sc_Me.new.ids)<-levels(my_sc_Me)
my_sc_Me<-RenameIdents(my_sc_Me,my_sc_Me.new.ids)

my_sc_Me.markers<-FindAllMarkers(my_sc_Me,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)
top10<-my_sc_Me.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
DefaultAssay(my_sc_Me)<-"RNA"
pdf(file = paste0(sc_path,"output/DoHeatmap_Me.pdf"))
DoHeatmap(my_sc_Me,features = top10$gene,size = 2,label=F)
dev.off()

pdf(file = paste0(sc_path,"output/UMAP/UMAPMe.pdf"))
DimPlot(my_sc_Me,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()


anno_Me<-data.frame(Barcode=colnames(my_sc_Me),Cluster=my_sc_Me@active.ident)
save(anno_Me,file = paste0(sc_path,"output/anno_Me.RData"))

################################output################################
load(paste0(sc_path,"output/anno_pre.RData"))
load(paste0(sc_path,"output/anno_m.RData"))
load(paste0(sc_path,"output/anno_e.RData"))
load(paste0(sc_path,"output/anno_h.RData"))
load(paste0(sc_path,"output/anno_Me.RData"))
anno_data<-rbind(anno_pre,anno_t,anno_m,anno_h,anno_Me,anno_e)

###filter clusters
anno_data$Cluster<-as.character(anno_data$Cluster)
anno_data<-subset(anno_data,anno_data$Cluster!="undefined")

anno_data$Cluster<-factor(anno_data$Cluster,levels = c("Hepatocyte","Endothelial","MDM","Mesenchymal","Cholangiocyte","T","DC",
                                                         "NK","Fibroblast","B","Kuffer","LEC","LSEC","Cyc"))

my_sc<-subset(my_sc,cells=as.vector(anno_data$Barcode))
anno_data<-anno_data[colnames(my_sc),]
Idents(my_sc)<-anno_data$Cluster
pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

DefaultAssay(my_sc)<-"RNA"
my_sc.markers<-FindAllMarkers(my_sc,only.pos = T)
my_sc.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_logFC)

top10<-my_sc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
write.table(top10,file = paste0(sc_path,"output/post_top10_markers.txt"),sep = "\t")

top5<-my_sc_f.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
DefaultAssay(my_sc)<-"RNA"
pdf(file = paste0(sc_path,"output/DoHeatmap.pdf"),height = 10,width = 15)
DoHeatmap(my_sc,features = top5$gene,size = 2,label = F)
dev.off()

save(anno_data,file = paste0(sc_path,"output/rolypoly/anno_data.RData"))
save(my_sc,file=paste0(sc_path,"output/my_sc.RData"))

################################block annotation################################
load(paste0(sc_path,"output/my_sc.RData"))
DefaultAssay(my_sc)<-"RNA"
library("biomaRt")
library(httr)
set_config(config(ssl_verifypeer = 0L))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",GRCh = 37)
ensembl_list<-read.table("/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE112271_multiple/13.a/genes.tsv"
                         ,header = F,stringsAsFactors = F)
gene_bp<-getBM(attributes = c("ensembl_gene_id","chromosome_name", "transcript_start","transcript_end","transcription_start_site","start_position","end_position","transcript_version"), 
               filters = "ensembl_gene_id", values = ensembl_list[,1], mart = ensembl)

c<-c(1:22)
c<-as.character(c)
gene_bp_filter<-gene_bp[which(gene_bp$chromosome_name %in% c),]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$ensembl_gene_id), ]


rownames(ensembl_list)<-ensembl_list[,1]
ensembl_list<-ensembl_list[gene_bp_filter$ensembl_gene_id,]
gene_bp_filter$external_gene_name<-ensembl_list[,2]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$external_gene_name), ]
###gene_coord
gene_coord_all<-gene_bp_filter[,c("external_gene_name","chromosome_name", "start_position","end_position")]
colnames(gene_coord_all)<-c("GENE","CHR","START","END")
rownames(gene_coord_all)<-gene_coord_all$GENE
write.table(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_coord/gene_coord_all.txt",
            sep = "\t",quote = F,col.names = T,row.names = F)

block_annotation<-gene_bp_filter
block_annotation$start_position<-block_annotation$transcription_start_site-10000
block_annotation$end_position<-block_annotation$transcription_start_site+10000
block_annotation<-block_annotation[, c("chromosome_name", "start_position","end_position","external_gene_name")]
colnames(block_annotation)<-c("chrom","start","end","label")
save(block_annotation,file =paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))


################################rolypoly input ################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/anno_data.RData"))

###gene_retain
counts_exp<-GetAssayData(my_sc,assay="RNA",slot="counts")
gene_retain<-CreateAssayObject(counts=counts_exp,min.cells = 3)%>%rownames()

###RNA matrix
RNA_exp<-GetAssayData(my_sc,assay="RNA",slot="data")[gene_retain,]%>%t()
RNA_exp<-as.data.frame(RNA_exp)
save(RNA_exp,file = paste0(sc_path,"output/rolypoly/RNA_exp.RData"))

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-anno_data
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined"
                                                      ,"Cyc","LSEC","Bud"))
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
                 cells=as.vector(anno_data[!anno_data$Cluster%in%c("undefined","Cyc","LSEC","Bud"),1]))
DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/10
DefaultAssay(my_sc_f)<-"RNA"
all_clusters<-c("Hepatocyte","Endothelial","MDM","Mesenchymal","Cholangiocyte","T","DC",
                "NK","Fibroblast","B","Kuffer","LEC")
###gene list
gene_set_findall<-FindAllMarkers(my_sc_f,
                                 min.pct=0.01,only.pos=T,return.thresh=1,logfc.threshold=0)

gene_set<-list()
for (c in 1:length(all_clusters)) {
  ident_x<-all_clusters[c]
  ident_y<-setdiff(all_clusters,ident_x)
  gene_set[[c]]<-FindMarkers(my_sc_f,ident.1=ident_x,
                             ident.2=ident_y,
                             logfc.threshold = 0,min.pct=0.01,
                             only.pos = T)
  gene_set_f<-gene_set[[c]]%>%top_n(n=-n_top,wt=p_val)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_list/EAS/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}



gene_coord_all<-read.table(file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_coord/gene_coord_all2.txt",
                           sep = "\t",header = T)
###conrol list
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_list/EAS/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)

###ldcts
path_list<-(paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/ldsc_out/EAS/",
                   1:length(all_clusters),
                   "_,/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/ldsc_out/EAS/0_"))

ldcts_data<-data.frame(all_clusters,path_list)
write.table(ldcts_data,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/GSE112271_multiple_EAS_ldcts"),
            sep = " ",quote = F,col.names =F,row.names =F)


