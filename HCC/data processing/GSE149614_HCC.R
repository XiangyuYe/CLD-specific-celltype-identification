################################SCT integrate################################
rm(list = ls())
gc()
library(scater)
library(Seurat)
library(dplyr)
library(scDblFinder)
library(bigreadr)
library(stringr)
options(future.globals.maxSize=4000*1024^2)
data_path="/net/mulan/disk2/yasheng/test/rolypoly/"
sc_path<-(paste0(data_path,"single_cell_data/GSE149614_HCC/"))

###data input
my_counts<-fread2(paste0(sc_path,"GSE149614_HCC_count.txt"))
rownames(my_counts)<-my_counts[,1]
my_counts<-my_counts[,-1]%>%as.matrix()
cell_list<-str_split(colnames(my_counts),"_")[,1]
tumor_list<-grep("T$",cell_list)
my_counts<-my_counts[,tumor_list]
patients<-cell_list[tumor_list]

###QC
pre_sc<-CreateSeuratObject(counts = my_counts)%>%as.SingleCellExperiment()
colData(pre_sc)<-DataFrame(patients)
pdf(paste0(sc_path,"output/1.pdf"))
pre_sc<-scDblFinder(pre_sc,verbose = T,clust.method = "overcluster")
dev.off()
table(pre_sc$scDblFinder.class)
my_sc<-pre_sc[,pre_sc$scDblFinder.class=="singlet"]%>%as.Seurat()

my_sc[["nCount_RNA"]] = colSums(x = my_sc, slot = "counts")  # nCount_RNA
my_sc[["nFeature_RNA"]] = colSums(x = GetAssayData(object = my_sc, slot = "counts") > 0)
my_sc[["percent.mt"]]<-PercentageFeatureSet(my_sc,pattern = "^MT-")

my_sc<-subset(my_sc,nFeature_RNA>=300&percent.mt<=20)
table(my_sc$patients)


###Normalize and Scale
my_sclist<-SplitObject(my_sc,split.by = "patients")
for (i in 1:length(my_sclist)) {
  my_sclist[[i]]<-SCTransform(my_sclist[[i]],
                              verbose = T, 
                              variable.features.n=3000
                              ,vars.to.regress = c("nCount_RNA","percent.mt"),
                              return.only.var.genes = F)
}

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

###PCA
DefaultAssay(my_sc)<-"integrated"
my_sc<-RunPCA(my_sc,npcs = 50)

pdf(file = paste0(sc_path,"output/DimPlot.pdf"))
DimPlot(my_sc,reduction = "pca",group.by  = "patients",ncol = 2)
dev.off()

tiff(file = paste0(sc_path,"output/DimHeatmap.tiff"),height = 5000,3000,compression = "lzw")
DimHeatmap(my_sc,dims = c(1:15),reduction = "pca",balanced = T)
dev.off()

pdf(file = paste0(sc_path,"output/Elbow.pdf"))
ElbowPlot(my_sc,ndims = 50)
dev.off()

my_sc<-JackStraw(my_sc,reduction = "pca",dims = 20)
my_sc<-ScoreJackStraw(my_sc,dims=1:20)
pdf(file = paste0(sc_path,"output/JackStraw.pdf"))
JackStrawPlot(my_sc,dims = 1:20)
dev.off()

###clustering
n_pca=40
resolution=1.2
my_sc<-FindNeighbors(my_sc,dims=1:n_pca)
my_sc<-FindClusters(my_sc,resolution = resolution)
my_sc<-RunUMAP(my_sc,dims = 1:n_pca)
 
pdf(file = paste0(sc_path,"output/UMAP/UMAP_g",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",group.by = "patients",label = T)
dev.off()

pdf(file = paste0(sc_path,"output/UMAP/UMAP",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",label = T)
dev.off()

table(my_sc@active.ident,my_sc$patients)

################################cell define################################
DefaultAssay(my_sc)<-"RNA"
my_sc<-NormalizeData(my_sc,verbose = T)
my_sc<-ScaleData(my_sc,vars.to.regress = c("nCount_RNA","percent.mt"))
save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))
#marker genes
T_features_list<-c("CD2","CD3D","CD3E","CD3G")
MP_features_list<-c("CD14","CD163","CD68","CSF1R")
B_features_list<-c("CD79A","MS4A1","JSRP1")
Endo_features_list<-c("ENG","CDH5","VWF")
Mesenchymal_features_list<-c("COL1A2","ACTA2","COL3A1","MYH11","RGS5")
Hep_features_list<-c("ALB","FGG","APOA1","APOC3")
Cyc_features_list<-c("MKI67","TOP2A","UBE2C","RRM2")
MAST_features_list<-c("TPSAB1","CPA3")
NK_features_list<-c("KLRF1","FGFBP2","TOX2","CD3D")

names_clusters<-c("T","MP","B","Endo","Mesenchymal","Hep","cyc","MAST","NK")
my_feature_list<-list(T_features_list,MP_features_list,B_features_list,Endo_features_list,Mesenchymal_features_list,Hep_features_list,Cyc_features_list,MAST_features_list,NK_features_list)

for (m in 1:length(names_clusters)) {
  pdf(file = paste0(sc_path,"output/markers/",names_clusters[[m]],"_post.pdf"))
  VlnPlot(my_sc,features =my_feature_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/FeaturePlot.pdf"))
FeaturePlot(my_sc,features =c("CD3D","ALB","CD68","ENG")
            ,pt.size = 0.001,label.size = 1)%>%print()
dev.off()

my_sc.markers<-FindMarkers(my_sc,ident.1=1,only.pos = T)
head(my_sc.markers,10)

my_sc.new.ids<-c("Hepatocyte",
                 "MP","Hepatocyte","Mixed","T","Hepatocyte","T","T","MP","Mesenchymal","T",
                 "Hepatocyte","MP","MP","Cyc","Plasma","Cyc","Hepatocyte","Hepatocyte","MP","MP",
                 "Hepatocyte","Cyc","Endothelial","B","Mast")
names(my_sc.new.ids)<-levels(my_sc)
my_sc<-RenameIdents(my_sc,my_sc.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_raw.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()

save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))

anno_pre<-data.frame(Barcode=colnames(my_sc),Cluster=my_sc@active.ident)
anno_pre<-subset(anno_pre,anno_pre$Cluster %in% c("Endothelial","B","Plasma","Mesenchymal","Cyc","Hepatocyte","Mast"))
anno_pre$Cluster<-factor(anno_pre$Cluster,ordered=F)
save(anno_pre,file = paste0(sc_path,"output/anno_pre.RData"))

###########Mixed########
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_x <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Mixed",1]))

###
DefaultAssay(my_sc_x)<-"integrated"
my_sc_x<-RunPCA(my_sc_x,npcs = 50)
pdf(file = paste0(sc_path,"output/Elbowx.pdf"))
ElbowPlot(my_sc_x,ndims = 50)
dev.off()

# my_sc_x<-JackStraw(my_sc_x,reduction = "pca",dims = 20)
# my_sc_x<-ScoreJackStraw(my_sc_x,dims=1:20)
# pdf(file = paste0(sc_path,"output/JackStraw_x.pdf"))
# JackStrawPlot(my_sc_x,dims = 1:20)
# dev.off()

###
n_pca=30
resolution=0.4
DefaultAssay(my_sc_x)<-"integrated"

my_sc_x<-FindNeighbors(my_sc_x,dims=1:n_pca)
my_sc_x<-FindClusters(my_sc_x,resolution = resolution)
my_sc_x<-RunUMAP(my_sc_x,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPx",n_pca,resolution,".pdf"))
DimPlot(my_sc_x,reduction = "umap",label = T)
dev.off()



DefaultAssay(my_sc_x)<-"RNA"

for (m in 1:length(names_clusters)) {
  pdf(file = paste0(sc_path,"output/markers/",names_clusters[[m]],"_x.pdf"))
  VlnPlot(my_sc_x,features =my_feature_list[[m]]
          ,pt.size = 0.01,ncol = 1)%>%print()
  dev.off()
}

my_sc.markers<-FindMarkers(my_sc,ident.1=2)
head(my_sc.markers,10)

pdf(file = paste0(sc_path,"output/markers/featurex.pdf"))
FeaturePlot(my_sc_x, c("ENG","VWF","APOC3","ALB"),pt.size=0.001,label.size=1)
dev.off()

my_sc_x.new.ids<-c("Hepatocyte",
                   "Hepatocyte","Endothelial")
names(my_sc_x.new.ids)<-levels(my_sc_x)
my_sc_x<-RenameIdents(my_sc_x,my_sc_x.new.ids)

anno_x<-data.frame(Barcode=colnames(my_sc_x),Cluster=my_sc_x@active.ident)
save(anno_x,file = paste0(sc_path,"output/anno_x.RData"))

###########T######
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
resolution=0.6
DefaultAssay(my_sc_t)<-"integrated"

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

names_T<-c("CD8","CD4","NK","cyc","B")
my_Tsub_list<-list(CD8_features_list,CD4_features_list,NK_features_list,Cyc_features_list,B_features_list)

for (m in 1:length(names_T)) {
  pdf(file = paste0(sc_path,"output/markers/",names_T[[m]],"_t.pdf"))
  VlnPlot(my_sc_t,features =my_Tsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

pdf(file = paste0(sc_path,"output/markers/featureT.pdf"))
FeaturePlot(my_sc_t, c("CD4", "CD8A", "KLRD1","CD79A"))
dev.off()

my_sc_t.markers<-FindMarkers(my_sc_t,min.pct = 0.25,ident.1=4)
head(my_sc_t.markers)

my_sc_t<-BuildClusterTree(my_sc_t,dims = 1:20)
pdf(file = paste0(sc_path,"output/BuildClusterTreeT.pdf"))
PlotClusterTree (object = my_sc_t)
dev.off()

my_sc_t.markers<-FindMarkers(my_sc_t,min.pct = 0.25,ident.1=1)
head(my_sc_t.markers)

###
my_sc_t.new.ids<-c("CD4+T",
                   "CD4+T","Treg","CD8+T","Plasma","Treg","NK","CD4+T","Cyc")
names(my_sc_t.new.ids)<-levels(my_sc_t)
my_sc_t<-RenameIdents(my_sc_t,my_sc_t.new.ids)

my_sc_t.markers<-FindAllMarkers(my_sc_t,only.pos = T)
top10_t<-my_sc_t.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap_T.pdf"))
DoHeatmap(my_sc_t,features = top10_t$gene,size = 2,label = F)
dev.off()

DefaultAssay(my_sc_t)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_t.pdf"))
DimPlot(my_sc_t,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

###
anno_t<-data.frame(Barcode=colnames(my_sc_t),Cluster=my_sc_t@active.ident)
save(anno_t,file=paste0(sc_path,"output/anno_t.RData"))

###########MP######
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_m <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="MP",1]))

###
DefaultAssay(my_sc_m)<-"integrated"
my_sc_m<-RunPCA(my_sc_m,npcs = 50)
pdf(file = paste0(sc_path,"output/Elbow_m.pdf"))
ElbowPlot(my_sc_m,ndims = 50)
dev.off()

# my_sc_m<-JackStraw(my_sc_m,reduction = "pca",dims = 20)
# my_sc_m<-ScoreJackStraw(my_sc_m,dims=1:20)
# pdf(file = paste0(sc_path,"output/JackStraw_m.pdf"))
# JackStrawPlot(my_sc_m,dims = 1:20)
# dev.off()

###
n_pca=20
resolution=0.6
DefaultAssay(my_sc_m)<-"integrated"

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

###
for (m in 1:length(names_MP)) {
  pdf(file = paste0(sc_path,"output/markers/",names_MP[[m]],"_m.pdf"))
  VlnPlot(my_sc_m,features =my_MPsub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc_m<-BuildClusterTree(my_sc_m,dims = 1:20)
pdf(file = paste0(sc_path,"output/BuildClusterTreeM.pdf"))
PlotClusterTree (object = my_sc_m)
dev.off()

my_sc.markers<-FindMarkers(my_sc_m,ident.1=5,only.pos=T)
head(my_sc.markers,10)

pdf(file = paste0(sc_path,"output/markers/FeaturePlotM.pdf"))
FeaturePlot(my_sc_m,features =c("C1QC","MNDA","TREM2","CD5L")
            ,pt.size = 0.001,label.size = 1)%>%print()
dev.off()

###
my_sc_m.new.ids<-c("Macrophage",
                   "Macrophage","Macrophage","Kuffer","DC","undefined","Macrophage","Macrophage","Macrophage","Kuffer","CLEC9A+DC")
names(my_sc_m.new.ids)<-levels(my_sc_m)
my_sc_m<-RenameIdents(my_sc_m,my_sc_m.new.ids)

my_sc_m.markers<-FindAllMarkers(my_sc_m,only.pos = T)
top10_m<-my_sc_m.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
DefaultAssay(my_sc_m)<-"RNA"
pdf(file = paste0(sc_path,"output/DoHeatmap_M.pdf"))
DoHeatmap(my_sc_m,features = top10_m$gene,size = 2,label=F)
dev.off()


DefaultAssay(my_sc_m)<-"integrated"
pdf(file = paste0(sc_path,"output/UMAP/UMAP_m.pdf"))
DimPlot(my_sc_m,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

###
anno_m<-data.frame(Barcode=colnames(my_sc_m),Cluster=my_sc_m@active.ident)
save(anno_m,file = paste0(sc_path,"output/anno_m.RData"))

###########Endo#####
ident_list <- data.frame(cell=names(Idents(my_sc)), cluster=Idents(my_sc))
my_sc_e <- subset(my_sc, cells=as.vector(ident_list[ident_list$cluster=="Mesenchymal",1]))


###
DefaultAssay(my_sc_e)<-"integrated"

my_sc_e<-RunPCA(my_sc_e,npcs = 50)
pdf(file = paste0(sc_path,"output/Elbow_e.pdf"))
ElbowPlot(my_sc_e,ndims = 50)
dev.off()

# my_sc_e<-JackStraw(my_sc_e,reduction = "pca",dims = 20)
# my_sc_e<-ScoreJackStraw(my_sc_e,dims=1:20)
# pdf(file = paste0(sc_path,"output/JackStraw_e.pdf"))
# JackStrawPlot(my_sc_e,dims = 1:20)
# dev.off()

###
n_pca=10
resolution=0.4

my_sc_e<-FindNeighbors(my_sc_e,dims=1:n_pca)
my_sc_e<-FindClusters(my_sc_e,resolution = resolution)
my_sc_e<-RunUMAP(my_sc_e,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAPe",n_pca,resolution,".pdf"))
DimPlot(my_sc_e,reduction = "umap",label = T)
dev.off()

###
DefaultAssay(my_sc_e)<-"RNA"

Fibro_features_list<-c("COL1A2","COL1A1")
SMC_features_list<-c("ACTA2")

names_e<-c("Fibro","SMC")
my_Esub_list<-list(Fibro_features_list,SMC_features_list)

for (m in 1:length(names_e)) {
  pdf(file = paste0(sc_path,"output/markers/",names_e[[m]],"_e.pdf"))
  VlnPlot(my_sc_e,features =my_Esub_list[[m]]
          ,pt.size = 0.1,ncol = 2)%>%print()
  dev.off()
}

my_sc_e<-BuildClusterTree(my_sc_e,dims = 1:10)
pdf(file = paste0(sc_path,"output/BuildClusterTreeE.pdf"))
PlotClusterTree (object = my_sc_e)
dev.off()

pdf(file = paste0(sc_path,"output/markers/FeaturePlotE.pdf"))
FeaturePlot(my_sc_e,features =c("COL1A2","COL1A1","ACTA2","ENG")
            ,pt.size = 0.001,label.size = 1)%>%print()
dev.off()


################################output################################
load(paste0(sc_path,"output/anno_m.RData"))
load(paste0(sc_path,"output/anno_t.RData"))
load(paste0(sc_path,"output/anno_x.RData"))
load(paste0(sc_path,"output/anno_pre.RData"))
anno_data<-rbind(anno_pre,anno_t,anno_m,anno_x)

anno_data$Cluster<-as.character(anno_data$Cluster)
anno_data<-subset(anno_data,anno_data$Cluster!="undefined")
anno_data$Cluster<-factor(anno_data$Cluster,levels = c("Hepatocyte","Macrophage","CD4+T","Mesenchymal","Treg","Monocyte",
                                                         "Plasma","CD8+T","DC","NK","Endothelial","B","Mast",
                                                         "Cyc"))


my_sc<-subset(my_sc,cells=as.vector(anno_data$Barcode))
anno_data<-anno_data[colnames(my_sc),]
Idents(my_sc)<-anno_data$Cluster

pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

DefaultAssay(my_sc)<-"RNA"
my_sc.markers<-FindAllMarkers(my_sc,only.pos = T)
top10<-my_sc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
write.table(top10,file = paste0(sc_path,"output/post_top10_markers.txt"),sep = "\t")

top5<-my_sc.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
pdf(file = paste0(sc_path,"output/DoHeatmap.pdf"),height = 10,width = 15)
DoHeatmap(my_sc,features = top5$gene,size = 2,label = F)
dev.off()

save(anno_data,file = paste0(sc_path,"output/rolypoly/anno_data.RData"))
save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))

################################block annotation#################################
load(paste0(sc_path,"output/my_sc.RData"))
ref_list<-read.table("/net/mulan/disk2/yasheng/test/rolypoly/single_cell_data/GSE136103_cirrhosis/cirrhosis/1_cd45+/genes.tsv"
                         ,header = F,stringsAsFactors = F)
rownames(ref_list)<-ref_list[,1]

DefaultAssay(my_sc)<-"RNA"
library(biomaRt)
library(httr)
ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",GRCh = 37)
ensembl_list<-subset(ref_list,ref_list[,2]%in%rownames(my_sc))
set_config(config(ssl_verifypeer = 0L))
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
save(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_coord_all.RData")
write.table(gene_coord_all,file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_coord/gene_coord_all.txt",
            sep = "\t",quote = F,col.names = T,row.names = F)


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
save(RNA_exp,file = paste0(sc_path,"output/rolypoly/RNA_exp.RData"))

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-anno_data
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc","Mast"))
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
anno_data2$Cluster[anno_data2$Cluster%in%c("Fibroblast")]<-"Mesenchymal"
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
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc","Mast"))
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
load(paste0(sc_path,"output/rolypoly/anno_TN.RData"))
anno_data<-anno_TN
my_sc_f<- subset(my_sc_f,
                 cells=as.vector(anno_data[!anno_data$Cluster%in%c("undefined","Cyc","Mast"),1]))

###gene_set
DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/10
DefaultAssay(my_sc)<-"RNA"
all_clusters<-c("Hepatocyte","Macrophage","CD4+T","Mesenchymal","Treg","Monocyte",
                "Plasma","CD8+T","DC","NK","Endothelial","B")

###gene list
gene_set_findall<-FindAllMarkers(my_sc_f,
                                 min.pct=0.001,only.pos=T,return.thresh=1,logfc.threshold=0)
save(gene_set_findall,file=paste0(sc_path,"output/rolypoly/gene_set_findall.RData"))

load(paste0(sc_path,"output/rolypoly/gene_set_findall.RData"))
for (c in 1:length(all_clusters)) {
  gene_set_f<-subset(gene_set_findall,gene_set_findall$cluster==all_clusters[c]&gene_set_findall$p_val_adj<0.05)
  gene_list<-rownames(gene_set_f)
  write.table(gene_list,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_list/sens/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}

###gene_ctrl
gene_coord_all<-read.table(file = "/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_coord/gene_coord_all.txt",
                           sep = "\t",header = T)
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/gene_list/sens/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)

###ldcts_data
path_list_FC<-(paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/ldsc_out/EAS/sens/",
                      1:length(all_clusters),
                      "_,/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/ldsc_out/EAS/sens/0_"))
ldcts_data<-data.frame(all_clusters,path_list_FC)
write.table(ldcts_data,file = paste0("/net/mulan/disk2/yasheng/test/LDSC_test/GSE149614_HCC/GSE149614_HCC_EAS_ldcts2"),
            sep = "\t",quote = F,col.names =F,row.names =F)

