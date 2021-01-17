################################data input################################
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
library(scDblFinder)
library(bigreadr)

options(future.globals.maxSize=4000*1024^2)
data_path="/home/integ_project/rolypoly/"
sc_path<-(paste0(data_path,"single_cell_data/GSE157783_IPDCO/"))

my_counts<-fread2(paste0(sc_path,"matrix.tsv"))
my_features<-fread2(paste0(sc_path,"features.tsv"))
my_barcodes<-fread2(paste0(sc_path,"barcodes.tsv"))

################################gene annotation################################
library("biomaRt")
library(httr)
set_config(config(ssl_verifypeer = 0L))

ensembl = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",GRCh = 37)
gene_list<-my_features$gene
gene_bp<-getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name", "transcript_start","transcript_end","transcription_start_site"), 
               filters = "ensembl_gene_id", values = gene_list, mart = ensembl)
c<-c(1:22)
c<-as.character(c)
gene_bp_filter<-gene_bp[which(gene_bp$chromosome_name %in% c),]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$ensembl_gene_id), ]
gene_bp_filter<-gene_bp_filter[!duplicated(gene_bp_filter$external_gene_name), ]

###gene_coord
gene_coord_all<-gene_bp_filter[,c("external_gene_name","chromosome_name", "transcript_start","transcript_end")]
colnames(gene_coord_all)<-c("GENE","CHR","START","END")
rownames(gene_coord_all)<-gene_coord_all$GENE
save(gene_coord_all,file = "/home/integ_project/LDSC_test/GSE157783_IPDCO/gene_coord_all.RData")
write.table(gene_coord_all,file = "/home/integ_project/LDSC_test/GSE157783_IPDCO/gene_coord/gene_coord_all.txt",
            sep = "\t",quote = F,col.names = T,row.names = F)

###block_annotation
block_annotation<-gene_bp_filter
block_annotation$start_position<-block_annotation$transcription_start_site-10000
block_annotation$end_position<-block_annotation$transcription_start_site+10000

block_annotation<-block_annotation[, c("chromosome_name", "start_position","end_position","external_gene_name")]
colnames(block_annotation)<-c("chrom","start","end","label")
save(block_annotation,file =paste0(sc_path,"output/rolypoly/block_annotation_all.RData"))

#######feature name transfer
rownames(my_features)<-my_features$gene
my_features<-my_features[,-1]
rownames(my_barcodes)<-my_barcodes$barcode
my_barcodes<-my_barcodes[,-1]
rownames(my_counts)<-rownames(my_features)
my_counts_f<-my_counts[gene_bp_filter$ensembl_gene_id,]
rownames(my_counts_f)<-gene_bp_filter$external_gene_name

pre_sc<-CreateSeuratObject(counts = my_counts_f)
pre_sc<-AddMetaData(pre_sc,my_barcodes)

pre_sc[["nCount_RNA"]] = colSums(x = pre_sc, slot = "counts")  # nCount_RNA
pre_sc[["nFeature_RNA"]] = colSums(x = GetAssayData(object = pre_sc, slot = "counts") > 0)
pre_sc[["percent.mt"]]<-PercentageFeatureSet(pre_sc,pattern = "^MT-")

my_sc<-subset(pre_sc,patient %in% paste0("PD",1:5))
my_sclist<-SplitObject(my_sc,split.by = "patient")

for (j in 1:5) {
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
                                   ,dims = 1:30,verbose = T)

my_sc<-IntegrateData(anchorset=my_anchors
                     ,normalization.method = "SCT"
                     ,dims = 1:30,verbose = T)


###Run PCA
DefaultAssay(my_sc)<-"integrated"
my_sc<-RunPCA(my_sc,npcs = 50)
###

pdf(file = paste0(sc_path,"output/Elbow.pdf"))
ElbowPlot(my_sc,ndims = 50)
dev.off()

###clustering
n_pca=50
resolution=0.8

my_sc<-FindNeighbors(my_sc,dims=1:n_pca)
my_sc<-FindClusters(my_sc,resolution = resolution)
my_sc<-RunUMAP(my_sc,dims = 1:n_pca)

pdf(file = paste0(sc_path,"output/UMAP/UMAP_g",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",group.by = "patient",label = T)
dev.off()
# 
pdf(file = paste0(sc_path,"output/UMAP/UMAP",n_pca,resolution,".pdf"))
DimPlot(my_sc,reduction = "umap",label = T)
dev.off()

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
Olig_features_list<-c("CNP","OLIG2","FA2H","MEGF11")
Miglia_features_list<-c("AIF1","CSF1R","CD53")
Astro_features_list<-c("GFAP","SOX9","ATP13A4","CBS")
Neuron_features_list<-c("RBFOX3","NRG1","KLHL1","CADPS2","GAD2","SLC17A6")
Endo_features_list<-c("VWF","ENG","CDH5","COL1A2","COL3A1","RGS5")
Peri_features_list<-c("RGS5","CSPG4","SLC6A12")
Epend_features_list<-c("DNAH11","DNAH9")
GABA_features_list<-c("GRIK1","SYNPR")
T_features_list<-c("CD2","KLRF1","CD3E","CD3G")

names_clusters<-c("Olig","Miglia","Astro","Neuron","Endo","Peri","Epend","GABA","T")
my_feature_list<-list(Olig_features_list,Miglia_features_list,Astro_features_list,Neuron_features_list,
                      Endo_features_list,Peri_features_list,Epend_features_list,GABA_features_list,T_features_list)

##volin plot
for (m in 1:length(names_clusters)) {
  pdf(file = paste0(sc_path,"output/markers/",names_clusters[[m]],"_post.pdf"))
  VlnPlot(my_sc,features =my_feature_list[[m]]
          ,pt.size = 0,ncol = 2)%>%print()
  dev.off()
}

##unsure cluster
my_sc.markers<-FindMarkers(my_sc,ident.1=16,only.pos = T)
my_sc.markers%>%top_n(10,wt=avg_logFC)

##cell type define
my_sc.new.ids<-c("Oligodendrocytes",
                 "Microglia","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Astrocytes","OPCs","Neuron","Astrocytes","Neuron","Endothelial",
                 "Pericytes","Ependymal","GABA","CADPS2+neurons","Fibroblast","T")
names(my_sc.new.ids)<-levels(my_sc)
my_sc<-RenameIdents(my_sc,my_sc.new.ids)

pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc,reduction = "umap",label = T,label.size =3 )+NoLegend()
dev.off()


################################seurat output################################
my_sc.markers<-FindAllMarkers(my_sc,only.pos = T)
top10<-my_sc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)
write.table(top10,file = paste0(sc_path,"output/top10_markers.txt"),sep = "\t")

top5<-my_sc.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC)
DefaultAssay(my_sc)<-"RNA"
pdf(file = paste0(sc_path,"output/DoHeatmap.pdf"),height = 10,width = 15)
DoHeatmap(my_sc,features = as.vector(top5$gene),size = 2,label = F)
dev.off()


pdf(file = paste0(sc_path,"output/UMAP/UMAP.pdf"))
DimPlot(my_sc_ff,reduction = "umap",label = T,label.size = 3)+NoLegend()
dev.off()

anno_data<-data.frame(Barcode=colnames(my_sc),Cluster=my_sc@active.ident)
save(anno_data,file = paste0(sc_path,"output/rolypoly/anno_data.RData"))
save(my_sc,file = paste0(sc_path,"output/my_sc.RData"))

################################rolypoly input################################
load(paste0(sc_path,"output/my_sc.RData"))
load(paste0(sc_path,"output/rolypoly/anno_data.RData"))
###gene_retain
counts_exp<-GetAssayData(my_sc,assay="RNA",slot="counts")
gene_retain<-CreateAssayObject(counts=counts_exp,min.cells = 3)%>%rownames()

###RNA matrix
RNA_exp<-GetAssayData(my_sc,assay="RNA",slot="data")%>%t()%>%as.data.frame()

gene_use<-intersect(as.character(gene_retain),
                    colnames(RNA_exp))
exp_data<-RNA_exp[,gene_use]
anno_data<-subset(anno_data,!anno_data$Cluster %in% c("undefined","Cyc","T","Fibroblast"))
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

DefaultAssay(my_sc_f)<-"SCT"
n_top<-nrow(my_sc_f)/10
DefaultAssay(my_sc_f)<-"RNA"
all_clusters<-c("Oligodendrocytes","Microglia","Astrocytes","Neuron","OPCs",
                "Endothelial","Pericytes","Ependymal","GABA","CADPS2+neurons")

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
              file = paste0("/home/integ_project/LDSC_test/GSE157783_IPDCO/gene_list/gene_list_",c,".txt"),
              sep = "\t",quote = F,col.names =F,row.names =F)
}

gene_coord_all<-read.table(file = "/home/integ_project/LDSC_test/GSE157783_IPDCO/gene_coord/gene_coord_all.txt",
                           sep = "\t",header = T)
###conrol list
gene_ctrl<-gene_coord_all$GENE
write.table(gene_ctrl,file = paste0("/home/integ_project/LDSC_test/GSE157783_IPDCO/gene_list/gene_list_0.txt"),
            sep = "\t",quote = F,col.names =F,row.names =F)

###ldcts
path_list<-(paste0("/home/integ_project/LDSC_test/GSE157783_IPDCO/ldsc_out/EAS/",
                   1:length(all_clusters),
                   "_,/home/integ_project/LDSC_test/GSE157783_IPDCO/ldsc_out/EAS/0_"))

ldcts_data<-data.frame(all_clusters,path_list)
write.table(ldcts_data,file = paste0("/home/integ_project/LDSC_test/GSE157783_IPDCO/GSE157783_IPDCO_EAS_ldcts"),
            sep = " ",quote = F,col.names =F,row.names =F)

      
      
