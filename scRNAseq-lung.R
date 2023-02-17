library(dplyr)
library(Matrix)
library(scran)
library(scater)
library(monocle3)
library(mpath)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(stringr)
library(AnnotationHub)
library(AnnotationDbi)
setwd("~\\lung")

dataan<-read.csv("lungexpressionmatrix.csv",row.names=1)

mc_dc <- CreateSeuratObject(counts =dataan, project = "pbmc2700", min.cells = 3, min.genes = 100)
dim(mc_dc@assays$RNA)
head(mc_dc@meta.data)
VlnPlot(mc_dc, features=c("nCount_RNA", "nFeature_RNA"), ncol = 2)  
FeatureScatter(object = mc_dc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
mc_dc@meta.data$countfeatureratio<-mc_dc@meta.data$nCount_RNA/mc_dc@meta.data$nFeature_RNA
mc_dc1<-subset(mc_dc,subset=nCount_RNA<4500000 & nFeature_RNA>1000)
mc_dc1<-subset(mc_dc1,subset=countfeatureratio<3000)
dim(mc_dc1@meta.data) 
VlnPlot(mc_dc1, features=c("nCount_RNA", "nFeature_RNA"), ncol = 2)  
FeatureScatter(object = mc_dc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
pbmc <- NormalizeData(mc_dc1,normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 2000)
pbmc1 <- ScaleData(pbmc1,vars.to.regress ="nCount_RNA")
pbmc2 <- RunPCA(pbmc1, features =VariableFeatures(object = pbmc1))
pbmc3 <- FindNeighbors(pbmc2, dims = 1:10)
pbmc4 <- FindClusters(pbmc3, resolution =0.7)
pbmc5 <- RunUMAP(pbmc4, dims = 1:10,check_duplicates = FALSE)
tsneplot<-UMAPPlot(pbmc5,label = TRUE,pt.size = 1.5)
tsneplot

##clustering of epithelial cells

pbmc5@meta.data$epithelial<-rep(0,dim(pbmc5@meta.data)[1])
pbmc5@meta.data$epithelial[which(pbmc5@meta.data$seurat_clusters==5)]<-rep(1,length(pbmc5@meta.data$epithelial[which(pbmc5@meta.data$seurat_clusters==5)]))
pbmc5@meta.data$epithelial[which(pbmc5@meta.data$seurat_clusters==6)]<-rep(1,length(pbmc5@meta.data$epithelial[which(pbmc5@meta.data$seurat_clusters==6)]))
pbmcci<-subset(pbmc5,subset=epithelial==1)
datai<-pbmcci@assays$RNA@counts
mc_dcc <- CreateSeuratObject(counts =datai, project = "pbmc2700", min.cells = 1, min.genes = 1)
pbmcc <- NormalizeData(mc_dcc,normalization.method = "LogNormalize", scale.factor = 10000)
pbmcc1 <- FindVariableFeatures(pbmcc,selection.method = "vst", nfeatures = 2000)
pbmcc1 <- ScaleData(pbmcc1,vars.to.regress ="nCount_RNA") 
pbmcc2 <- RunPCA(pbmcc1, features =VariableFeatures(object = pbmcc1))
pbmcc3 <- FindNeighbors(pbmcc2, dims = 1:10)
pbmcc4 <- FindClusters(pbmcc3, resolution =0.7)
pbmcc5 <- RunUMAP(pbmcc4, dims = 1:10,check_duplicates = FALSE)
tsneplott<-UMAPPlot(pbmcc5,label = TRUE,pt.size = 1.5)
tsneplott
pbmcc5$AGER<-pbmcc5@assays$RNA@data["AGER",]
pbmcc5$SFTPC<-pbmcc5@assays$RNA@data["SFTPC",]
pbmcc5$TPPP3<-pbmcc5@assays$RNA@data["TPPP3",]
pbmcc5$SCGB3A2<-pbmcc5@assays$RNA@data["SCGB3A2",]
pbmcc5@meta.data$celltype<-rep("Other epithelial",dim(pbmcc5@meta.data)[1])
AT1=intersect(which(pbmcc5@meta.data$AGER>1),which(exp(pbmcc5@meta.data$AGER)/length["AGER",5] > exp(pbmcc5@meta.data$SFTPC)/length["SFTPC",5]))
pbmcc5@meta.data$celltype[AT1]<-rep("AT1",length(AT1))
AT2=intersect(which(pbmcc5@meta.data$SFTPC>1),which(exp(pbmcc5@meta.data$SFTPC)/length["SFTPC",5] > exp(pbmcc5@meta.data$AGER)/length["AGER",5]))
pbmcc5@meta.data$celltype[AT2]<-rep("AT2",length(AT2))
CT=which(pbmcc5@meta.data$TPPP3>3)
pbmcc5@meta.data$celltype[CT]<-rep("Ciliated cells",length(CT))

#cluster of T cells

pbmc5@meta.data$tcell<-rep(0,dim(pbmc5@meta.data)[1])
pbmc5@meta.data$tcell[which(pbmc5@meta.data$seurat_clusters==2)]<-rep(1,length(pbmc5@meta.data$seurat_clusters[which(pbmc5@meta.data$seurat_clusters==2)]))
pbmce<-subset(pbmc5,subset=tcell==1)
datae<-pbmce@assays$RNA@counts
mc_de <- CreateSeuratObject(counts =datae, project = "pbmc2700", min.cells = 1, min.genes = 1)
pbmce <- NormalizeData(mc_de,normalization.method = "LogNormalize", scale.factor = 10000)
pbmce1 <- FindVariableFeatures(pbmce,selection.method = "vst", nfeatures = 2000)
pbmce1 <- ScaleData(pbmce1,vars.to.regress ="nCount_RNA") 
pbmce2 <- RunPCA(pbmce1, features =VariableFeatures(object = pbmce1))
pbmce3 <- FindNeighbors(pbmce2, dims = 1:10)
pbmce4 <- FindClusters(pbmce3, resolution =0.7)
pbmce5 <- RunUMAP(pbmce4, dims = 1:10,check_duplicates = FALSE)
tsneplote<-UMAPPlot(pbmce5,label = TRUE,pt.size = 1.5)
tsneplote

pbmce5$CD4<-pbmce5@assays$RNA@data["CD4",]
pbmce5$CD8A<-pbmce5@assays$RNA@data["CD8A",]
pbmce5$CD8B<-pbmce5@assays$RNA@data["CD8B",]
pbmce5@meta.data$celltype<-rep("Other T",dim(pbmce5@meta.data)[1])

CD8=union(which(pbmce5@meta.data$CD8A>0.5),which(pbmce5@meta.data$CD8B>0.5))
CD4=which(pbmce5@meta.data$CD4>0.5)
pbmce5@meta.data$celltype[setdiff(CD8,CD4)]<-rep("CD8+T",length(setdiff(CD8,CD4)))
pbmce5@meta.data$celltype[setdiff(CD4,CD8)]<-rep("CD4+T",length(setdiff(CD4,CD8)))

##cluster of monocyte-macrophage

pbmc5@meta.data$macromono<-rep(0,dim(pbmc5@meta.data)[1])
pbmc5@meta.data$macromono[which(pbmc5@meta.data$seurat_clusters==1)]<-rep(1,length(pbmc5@meta.data$seurat_clusters[which(pbmc5@meta.data$seurat_clusters==1)]))
pbmc5@meta.data$macromono[which(pbmc5@meta.data$seurat_clusters==4)]<-rep(1,length(pbmc5@meta.data$seurat_clusters[which(pbmc5@meta.data$seurat_clusters==4)]))
pbmcmi<-subset(pbmc5,subset=macromono==1) 
datam<-pbmcmi@assays$RNA@counts
mc_dcm <- CreateSeuratObject(counts =datam, project = "pbmc2700", min.cells = 1, min.genes = 1)
pbmcm <- NormalizeData(mc_dcm,normalization.method = "LogNormalize", scale.factor = 10000)
pbmcm1 <- FindVariableFeatures(pbmcm,selection.method = "vst", nfeatures = 2000)
pbmcm1 <- ScaleData(pbmcm1,vars.to.regress ="nCount_RNA")
pbmcm2 <- RunPCA(pbmcm1, features =VariableFeatures(object = pbmcm1))
pbmcm3 <- FindNeighbors(pbmcm2, dims = 1:20)
pbmcm4 <- FindClusters(pbmcm3, resolution =0.8)
pbmcm5 <- RunUMAP(pbmcm4, dims = 1:20,check_duplicates = FALSE)
tsneplott<-UMAPPlot(pbmcm5,label = TRUE,pt.size = 1.5)
tsneplott
pbmcm5$FCGR3B<-pbmcm5@assays$RNA@data["FCGR3B",]
pbmcm5$FCN1<-pbmcm5@assays$RNA@data["FCN1",]
pbmcm5@meta.data$celltype<-rep("macrophage",dim(pbmcm5@meta.data)[1])
DEN=which(pbmcm5@meta.data$seurat_clusters==1)
pbmcm5@meta.data$celltype[DEN]<-rep("Dendritic cells",length(DEN))
MONO=union(which(pbmcm5@meta.data$seurat_clusters==0),intersect(which(pbmcm5@meta.data$seurat_clusters==4),which(pbmcm5@meta.data$FCN1>1)))
pbmcm5@meta.data$celltype[MONO]<-rep("Monocytes",length(MONO))
GRA=which(pbmcm5@meta.data$FCGR3B>0.8)
pbmcm5@meta.data$celltype[GRA]<-rep("Granulocyte",length(GRA))



new.cluster.ids<-c("Endothelial","Myeloid","T cell","Fibroblast","Myeloid","Epithelial","Epithelial","Mast cell","Fibroblast","Fibroblast","NK cell","B cell")
names(new.cluster.ids) <- levels(pbmc5)
pbmc6 <- RenameIdents(pbmc5,new.cluster.ids)
tsneplot2<-UMAPPlot(pbmc6,label = TRUE,pt.size = 2)
tsneplot2
pbmc6@meta.data$celltype1<-as.character(pbmc6@active.ident)
pbmc6@meta.data$celltype2<-as.character(pbmc6@active.ident)

for (i in which(pbmc6@meta.data$epithelial==1))
{
  pbmc6@meta.data$celltype2[i]<-pbmcc5@meta.data$celltype[which(rownames(pbmcc5@meta.data)==rownames(pbmc6@meta.data)[i])]
}
for (i in which(pbmc6@meta.data$tcell==1))
{
  pbmc6@meta.data$celltype2[i]<-pbmce5@meta.data$celltype[which(rownames(pbmce5@meta.data)==rownames(pbmc6@meta.data)[i])]
}

for (i in which(pbmc6@meta.data$macromono==1))
{
  pbmc6@meta.data$celltype2[i]<-pbmcm5@meta.data$celltype[which(rownames(pbmcm5@meta.data)==rownames(pbmc6@meta.data)[i])]
}

pbmc6@meta.data$celltype1[!duplicated(pbmc6@meta.data$celltype1)]
pbmc6@meta.data$celltype2[!duplicated(pbmc6@meta.data$celltype2)]

lungcell<-pbmc6@meta.data[,c(12,13)]
setwd("~\\lung")
write.csv(lungcell,"lungcellgroup.csv")





setwd("~\\lung")
group<-read.csv("lungcellgroup.csv")
identical(group[,1],rownames(pbmc6@meta.data))
pbmc6$celltype2<-group$celltype2
pbmc6@active.ident<-as.factor(pbmc6$celltype2)
names(pbmc6@active.ident)<-rownames(pbmc6@meta.data)
UMAPPlot(pbmc6,pt.size=1.5)

#markers of each group：
marker2gene<-c("AGER","SFTPC","CD19","CD4","CD8A","TPPP3","CD1C","PECAM1","PDGFRA","CSF3R","CSF1R","CPA3","CYBB","KLRD1","EPCAM","CD3D")
marker2cell<-c("AT1","AT2","B cell","CD4+T","CD8+T","Ciliated cells","Dendritic cells","Endothelial","Fibroblast","Granulocyte","macrophage","Mast cell","Monocytes","NK cell","Other epithelial","Other T")

 library(reshape2)
vln.df=as.data.frame(pbmc6[["RNA"]]@data[marker2gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
pbmc6@meta.data$CB<-rownames(pbmc6@meta.data)
anno=pbmc6@meta.data[,c("CB","celltype2")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = marker2gene) 
vln.df$celltype=factor(vln.df$celltype,levels = marker2cell)
ggplot(vln.df,aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(gene~.,scales = "free_y")+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 180,hjust = 1,vjust = 1),strip.text=element_text(size=35),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none",axis.text.x=element_text(size=50),axis.text.y=element_text(size=25),axis.title.y=element_text(size=50)
  )
ggsave("violinplot.png",width = 40,height = 65,limitsize=FALSE)

pbmc6<-pbmc5
setwd("~\\lung")
group<-read.csv("lungcellgroup.csv")
identical(group[,1],rownames(pbmc6@meta.data))
pbmc6$celltype2<-group$celltype2
marker2cell<-c("AT1","AT2","B cell","CD4+T","CD8+T","Ciliated cells","Dendritic cells","Endothelial","Fibroblast","Granulocyte","macrophage","Mast cell","Monocytes","NK cell","Other epithelial","Other T")
numbercell<-c()
for (i in c(1:length(marker2cell)))
{
  numbercell[i]<-length(which(pbmc6@meta.data$celltype2==marker2cell[i]))
}
names(numbercell)<-marker2cell
numberc<-data.frame(marker2cell,numbercell)
numberc$ratio<-numberc$number/sum(numberc$number)
colnames(numberc)<-c("celltype","number","ratio")
gg<-ggplot(numberc,aes(x=celltype,y=number))+geom_bar(aes(fill=celltype),stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle = 45,size=25,vjust=0.5),axis.title.x=element_text(size=35),axis.text.y=element_text(size=35),axis.title.y=element_text(size=25),legend.position = "none")
ggsave("cellnumber.png",gg,width = 40,height = 12,limitsize=FALSE)
gg<-ggplot(numberc,aes(x=celltype,y=ratio))+geom_bar(aes(fill=celltype),stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle = 45,size=25,vjust=0.5),axis.title.x=element_text(size=35),axis.text.y=element_text(size=35),axis.title.y=element_text(size=25),legend.position = "none")
ggsave("cellratio.png",gg,width = 40,height = 12,limitsize=FALSE)
