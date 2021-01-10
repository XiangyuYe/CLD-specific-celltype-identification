################################SCTransform integrate################################
rm(list = ls())
gc()
library(scater)
library(Seurat)
library(dplyr)
library(scDblFinder)
library(bigreadr)
options(future.globals.maxSize=4000*1024^2)
data_path="/net/mulan/disk2/yasheng/test/rolypoly/"
sc_path<-"/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/integrated_HCC/"
sc_path1<-"/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE149614_HCC/"
sc_path2<-"/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE112271_multiple/"

###TN data
load(paste0(sc_path1,"output/my_sc.RData"))
my_sclist1<-SplitObject(my_sc,split.by = "patients")
for (i in 1:length(my_sclist)) {
  my_sclist1[[i]]<-SCTransform(my_sclist[[i]],
                               verbose = T, 
                               variable.features.n=3000
                              ,vars.to.regress = c("nCount_RNA","percent.mt"),
                              return.only.var.genes = F
  )
}

###multile data
load(paste0(sc_path2,"output/my_sc.RData"))
my_sclist2<-SplitObject(my_sc,split.by = "orig.ident")

for (i in 1:length(my_sclist)) {
  my_sclist2[[i]]<-SCTransform(my_sclist[[i]],
                               verbose = T, 
                               variable.features.n=3000,
                               vars.to.regress = c("nCount_RNA","percent.mt"),
                               return.only.var.genes = F)
}

###integrate
my_sclist<-c(my_sclist1,my_sclist2)
my.features<-SelectIntegrationFeatures(object.list =my_sclist
                                       ,nfeatures = 3000)

my_sclist<-PrepSCTIntegration(object.list = my_sclist
                              ,anchor.features = my.features
                              ,verbose = T)

my_anchors<-FindIntegrationAnchors(object.list=my_sclist
                                   ,normalization.method = "SCT"
                                   ,anchor.features = my.features
                                   ,dims = 1:30,verbose = T)

my_sc<-IntegrateData(anchorset=my_anchors
                     ,normalization.method = "SCT"
                     ,dims = 1:30,verbose = T)

###Run PCA
my_sc<-RunPCA(my_sc,npcs = 80)

percent.var <- my_sc@ reductions$ pca@ stdev
print(paste0("sum_var30 is ",sum(percent.var[1:30])))

pdf(file = paste0(sc_path,"output/DimPlot.pdf"))
DimPlot(my_sc,reduction = "pca",group.by  = "orig.ident")%>%print()
dev.off()

pdf(file = paste0(sc_path,"output/Elbow.pdf"))
ElbowPlot(my_sc,ndims = 80)%>%print()
dev.off()

###clustering
n_pca=80
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
my_sc.markers<-FindMarkers(my_sc,ident.1=22,ident.2=5,only.pos = T)
head(my_sc.markers,10)

##cell type define
my_sc.new.ids<-c("Epithelial",
                 "Epithelial","MP","Epithelial","MP","Mesenchyme","T","Endothelial","T","Epithelial","Epithelial",
                 "MP","Epithelial","Epithelial","Endothelial","Epithelial","Cyc","Cyc","Epithelial","Plasma","B",
                 "Endothelial","MP")
names(my_sc.new.ids)<-levels(my_sc)
my_sc<-RenameIdents(my_sc,my_sc.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_raw.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_pre<-data.frame(Barcode=colnames(my_sc),Cluster=my_sc@active.ident)
anno_pre<-subset(anno_pre,anno_pre$Cluster %in% c("Epithelial","B","Plasma","Cyc","Endothelial"))
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
n_pca=30
resolution=0.8
my_sc_t<-FindNeighbors(my_sc_t,dims=1:n_pca)
my_sc_t<-FindClusters(my_sc_t,resolution = resolution)
my_sc_t<-RunUMAP(my_sc_t,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPt",n_pca,resolution,".pdf"))
DimPlot(my_sc_t,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_t)<-"RNA"

CD8_features_list<-c("CD8A", "CD8B","GZMA","GZMH")
CD4_features_list<-c("CD4","CCR6","CCR7","SELL","FOXP3","RORC","IL17A")
NK_features_list<-c("KLRF1","FGFBP2","TOX2","CD3D")

names_T<-c("CD8","CD4","NK","cyc","Mast")
my_Tsub_list<-list(CD8_features_list,CD4_features_list,NK_features_list,Cyc_features_list,MAST_features_list)

for (m in 1:length(names_T)) {
  pdf(file = paste0(sc_path,"output/markers/",names_T[[m]],"_t.pdf"))
  VlnPlot(my_sc_t,features =my_Tsub_list[[m]]
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/featureT.pdf"))
FeaturePlot(my_sc_t, c("CD4", "CD8A", "KLRF1","FOXP3"))
dev.off()

my_sc.markers<-FindMarkers(my_sc_t,ident.1=11,only.pos = T)
head(my_sc.markers,10)

my_sc_t.new.ids<-c("CD4+T",
                   "CD4+T","NK","CD8+T","NK","CD8+T","CD4+T","CD4+T","Treg","Treg","CD4+T",
                   "Treg")
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
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/featureM.pdf"))
FeaturePlot(my_sc_m, c("S100A8", "FCN1", "CD1E","C1QC"))
dev.off()

my_sc.markers<-FindMarkers(my_sc_m,ident.1=14,only.pos = T)
head(my_sc.markers,10)

my_sc_m.new.ids<-c("Macrophage",
                   "Macrophage","Macrophage","Macrophage","Monocyte","DC","Macrophage","Macrophage","Monocyte","Macrophage","Macrophage",
                   "Macrophage","Macrophage","Macrophage","undefined","Macrophage")
names(my_sc_m.new.ids)<-levels(my_sc_m)
my_sc_m<-RenameIdents(my_sc_m,my_sc_m.new.ids)

DefaultAssay(my_sc_m)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_m.pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_m<-data.frame(Barcode=colnames(my_sc_m),Cluster=my_sc_m@active.ident)
save(anno_m,file = paste0(sc_path,"output/anno_m.RData"))
###########Endo cells###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_e <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Endothelial",1]))
DefaultAssay(my_sc_e)

###
DefaultAssay(my_sc_e)<-"integrated"
my_sc_e<-RunPCA(my_sc_e,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_e.pdf"))
ElbowPlot(my_sc_e,ndims = 50)
dev.off()

###
n_pca=20
resolution=0.4
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
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

my_sc.markers<-FindMarkers(my_sc_e,ident.1=7,only.pos = T)
head(my_sc.markers,10)

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

# my_sc_h<-JackStraw(my_sc_h,reduction = "pca",dims = 20)
# my_sc_h<-ScoreJackStraw(my_sc_h,dims=1:20)
# pdf(file = paste0(sc_path,"output/JackStraw_h.pdf"))
# JackStrawPlot(my_sc_h,dims = 1:20)
# dev.off()

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



###########Mesenchymal cells subdivide###########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_p <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Mesenchyme",1]))
DefaultAssay(my_sc_p)
###
DefaultAssay(my_sc_p)<-"integrated"
my_sc_p<-RunPCA(my_sc_p,npcs = 50)

pdf(file = paste0(sc_path,"output/Elbow_p.pdf"))
ElbowPlot(my_sc_p,ndims = 50)
dev.off()

###
n_pca=20
resolution=0.4

my_sc_p<-FindNeighbors(my_sc_p,dims=1:n_pca)
my_sc_p<-FindClusters(my_sc_p,resolution = resolution)
my_sc_p<-RunUMAP(my_sc_p,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPp",n_pca,resolution,".pdf"))
DimPlot(my_sc_p,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_p)<-"RNA"
pdf(file = paste0(sc_path,"output/markers/featureH.pdf"))
FeaturePlot(my_sc_p, Mesenchymal_features_list[1:4])
dev.off()

my_sc_p.new.ids<-c("Mesenchyme",
                   "Mesenchyme","Mesenchyme","Mesenchyme","Fibroblast","Mesenchyme")
names(my_sc_p.new.ids)<-levels(my_sc_p)
my_sc_p<-RenameIdents(my_sc_p,my_sc_p.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAPp.pdf"))
DimPlot(my_sc_p,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

anno_p<-data.frame(Barcode=colnames(my_sc_p),Cluster=my_sc_p@active.ident)
save(anno_p,file = paste0(sc_path,"output/anno_p.RData"))

################################output################################
load(paste0(sc_path,"output/anno_pre.RData"))
load(paste0(sc_path,"output/anno_t.RData"))
load(paste0(sc_path,"output/anno_m.RData"))
load(paste0(sc_path,"output/anno_h.RData"))
load(paste0(sc_path,"output/anno_p.RData"))
load(paste0(sc_path,"output/my_sc.RData"))
anno_data<-rbind(anno_pre,anno_t,anno_m,anno_p,anno_h)
###filter and order
anno_data<-subset(anno_data,anno_data$Cluster!="undefined")

anno_data$Cluster<-factor(anno_data$Cluster,levels = c("Epithelial","Macrophage","Endothelial","CD4+T","Mesenchyme","NK",
                                                       "Monocyte","CD8+T","DC","Treg","Plasma","Fibroblast",
                                                       "B","Cyc"))
my_sc<-subset(my_sc,cells=as.vector(anno_data$Barcode))
anno_data<-anno_data[colnames(my_sc),]
Idents(my_sc)<-anno_data$Cluster

pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))
###marker gene present
DefaultAssay(my_sc)<-"RNA"
my_sc.markers<-FindAllMarkers(my_sc,only.pos = T)
my_sc.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_logFC)

top10<-my_sc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
write.table(top10,file = paste0(sc_path,"output/post_top10_markers_f.txt"),sep = "\t")

top5<-my_sc.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap.pdf"),height = 10,width = 15)
DoHeatmap(my_sc,features = top5$gene,size = 2,label = F)
dev.off()

###data output
save(anno_data,file = paste0(sc_path,"output/rolypoly/anno_data.RData"))
save(my_sc,file = paste0(sc_path,"output/rolypoly/my_sc.RData"))
################################block annotation#################################
load(paste0(sc_path,"output/my_sc.RData"))
gene_coord1<-read.table("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_coord/gene_coord_all.txt",header = T)
gene_coord2<-read.table("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_coord/gene_coord_all.txt",header = T)
gene_coord_all<-rbind(gene_coord1,gene_coord2)
gene_coord_all<-gene_coord_all[!duplicated(gene_coord_all$GENE), ]

save(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_coord_all.RData")
write.table(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_coord/gene_coord_all.txt",
            sep = "\t",quote = F,col.names = T,row.names = F)

load(paste0(data_path,"single_cell_data/GSE149614_HCC/output/rolypoly/block_annotation_all.RData"))
block_annotation1<-block_annotation
load(paste0(data_path,"single_cell_data/GSE112271_multiple/output/rolypoly/block_annotation_all.RData"))
block_annotation2<-block_annotation
block_annotation<-rbind(block_annotation1,block_annotation2)
block_annotation<-block_annotation[!duplicated(block_annotation$label), ]

save(block_annotation,file =paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))



################################rolypoly input################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/anno_data.RData"))
load(paste0(sc_path,"output/rolypoly/RNA_exp.RData"))
load(paste0(sc_path,"output/rolypoly/gene_retain.RData"))
###gene_retain
counts_exp<-GetAssayData(my_sc,assay="RNA",slot="counts")
gene_retain<-CreateAssayObject(counts=counts_exp,min.cells = 3)%>%rownames()
save(gene_retain,file=paste0(sc_path,"output/rolypoly/gene_retain.RData"))
###RNA matrix
RNA_exp<-GetAssayData(my_sc,assay="RNA",slot="data")[gene_retain,]%>%t()
RNA_exp<-as.data.frame(RNA_exp)
save(RNA_exp,file = paste0(sc_path,"output/rolypoly/RNA_exp.RData"))

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-anno_data
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc"))
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

###sensitive analysis
anno_data2<-anno_data
anno_data2$Cluster<-as.character(anno_data2$Cluster)
anno_data2$Cluster[anno_data2$Cluster%in%c("MDM","Kuffer")]<-"Macrophage"
anno_data2$Cluster[anno_data2$Cluster%in%c("LEC","LSEC")]<-"Endothelial"
anno_data2$Cluster[anno_data2$Cluster%in%c("Fibroblast","Mesenchyme")]<-"Mesenchymal"
anno_data2$Cluster[anno_data2$Cluster%in%c("CD4+T","Treg","CD8+T")]<-"T"
anno_data2$Cluster[anno_data2$Cluster%in%c("Plasma")]<-"B"
anno_data2$Cluster[anno_data2$Cluster%in%c("NKT","cNK")]<-"NK"
anno_data2$Cluster[anno_data2$Cluster%in%c("pDC","CLEC9A+DC")]<-"DC"
anno_data2$Cluster[anno_data2$Cluster%in%c("Hepatocyte","Cholangiocyte")]<-"Epithelial"

save(anno_data2,file = paste0(sc_path,"output/rolypoly/anno_data2.RData"))

gene_retain<-CreateAssayObject(counts=counts_exp,min.cells = 3)%>%rownames()
load(paste0(sc_path,"output/rolypoly/RNA_exp.RData"))

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-anno_data2
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc","LSEC","Bud"))
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

cell_exp2<-abs(exp_data)
save(cell_exp2,file = paste0(sc_path,"output/rolypoly/cell_exp2.RData"))

################################LDSC_cts input################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))
my_sc_f<-my_sc[block_annotation$label,]

# load(paste0(sc_path,"output/anno_data.RData"))
my_sc_f<- subset(my_sc_f,
                 cells=as.vector(anno_data[!anno_data$Cluster%in%c("undefined","Cyc"),1]))
DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/10
DefaultAssay(my_sc_f)<-"RNA"

###gene list
all_clusters<-unique(anno_data$Cluster)

###gene list
gene_set_findall<-FindAllMarkers(my_sc_f,
                                 min.pct=0.001,only.pos=T,return.thresh=1,logfc.threshold=0)
save(gene_set_findall,paste0(sc_path,"output/rolypoly/gene_set_findall.RData"))

for (c in 1:length(all_clusters)) {
  gene_set_f<-subset(gene_set_findall,gene_set_findall$cluster==all_clusters[c])
  gene_set_f<-gene_set[[c]]%>%top_n(n=-n_top,wt=p_val)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE112271_multiple/gene_list/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}

###gene_ctrl
gene_coord_all<-read.table(file = "/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_coord/gene_coord_all.txt",
                           sep = "\t",header = T)
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_list/EAS/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)


###ldcts
path_list<-(paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/ldsc_out/EAS/",
                   1:length(all_clusters),
                   "_,/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/ldsc_out/EAS/0_"))

ldcts_data<-data.frame(all_clusters,path_list)
write.table(ldcts_data,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/integrated_HCC_EAS_ldcts"),
            sep = " ",quote = F,col.names =F,row.names =F)


###sensitive analysis
load(paste0(sc_path,"output/rolypoly/gene_set_findall.RData"))
all_clusters<-c("Epithelial","Macrophage","Endothelial","CD4+T","Mesenchyme","NK",
                "Monocyte","CD8+T","DC","Plasma","Fibroblast",
                "B","Treg")
n_top<-nrow(my_sc_f)/20
for (c in 1:length(all_clusters)) {
  gene_set_f<-subset(gene_set_findall,gene_set_findall$cluster==all_clusters[c])
  gene_set_f<-gene_set[[c]]%>%top_n(n=-n_top,wt=p_val)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_list/sens/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}

###conrol list
gene_coord_all<-read.table(file = "/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_coord/gene_coord_all2.txt",
                           sep = "\t",header = T)
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/gene_list/sens/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)
###ldcts
path_list<-(paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/ldsc_out/EAS/sens/",
                   1:length(all_clusters),
                   "_,/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/ldsc_out/EAS/sens/0_"))

ldcts_data<-data.frame(all_clusters,path_list)
write.table(ldcts_data,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/integrated_HCC/sens/integrated_HCC_EAS_ldcts"),
            sep = " ",quote = F,col.names =F,row.names =F)

#####
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/anno_data2.RData"))
load(paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))
anno_data2<-anno_data2[colnames(my_sc),]
Idents(my_sc)<-anno_data2$Cluster
my_sc_f<-my_sc[block_annotation$label,]

# load(paste0(sc_path,"output/anno_data.RData"))
my_sc_f<- subset(my_sc_f,
                 cells=as.vector(anno_data2[!anno_data2$Cluster%in%c("undefined","Cyc"),1]))
DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/20
DefaultAssay(my_sc_f)<-"RNA"

all_clusters<-unique(anno_data2$Cluster)

gene_set_findall2<-FindAllMarkers(my_sc_f,
                                 min.pct=0.001,only.pos=T,return.thresh=1,logfc.threshold=0)
save(gene_set_findall2,file=paste0(sc_path,"output/rolypoly/gene_set_findall2.RData"))

load(paste0(sc_path,"output/rolypoly/gene_set_findall2.RData"))
for (c in 1:length(all_clusters)) {
  gene_set_f<-subset(gene_set_findall2,gene_set_findall2$cluster==all_clusters[c]&gene_set_findall2$p_val_adj<0.05)
  gene_set_f<-gene_set[[c]]%>%top_n(n=-n_top,wt=p_val)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_list/sens/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}
