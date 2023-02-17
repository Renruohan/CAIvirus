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
setwd("~\\PBMC")

dataan<-read.csv("PBMCexpressionmatrix.csv",row.names=1)
dataan<-dataan[,-c(1:5)]

mc_dc <- CreateSeuratObject(counts =dataan, project = "pbmc2700", min.cells = 3, min.genes = 200)
dim(mc_dc@assays$RNA)
mc_dc[["percent.mt"]] <- PercentageFeatureSet(mc_dc, pattern = "^MT-")
head(mc_dc@meta.data)
VlnPlot(mc_dc, features=c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)  

par(mfrow = c(1, 2))
FeatureScatter(object = mc_dc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = mc_dc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

mc_dc1<-subset(mc_dc,subset=nCount_RNA>1000 & nFeature_RNA>100 & percent.mt<4) #13740
dim(mc_dc1@meta.data)
VlnPlot(mc_dc1, features=c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)  
pbmc <- NormalizeData(mc_dc1,normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 1000)
pbmc1 <- ScaleData(pbmc1,vars.to.regress ="nCount_RNA")
pbmc2 <- RunPCA(pbmc1, features =VariableFeatures(object = pbmc1))
pbmc3 <- FindNeighbors(pbmc2, dims = 1:50)
pbmc4 <- FindClusters(pbmc3, resolution =1)
pbmc5 <- RunUMAP(pbmc4, dims = 1:50,check_duplicates = FALSE)
tsneplot<-UMAPPlot(pbmc5,label = TRUE,pt.size = 1.5)
tsneplot
pbmc5@meta.data$seurat_clusters<-as.character(pbmc5@meta.data$seurat_clusters)

l<-union(which(pbmc5@assays$RNA@data["ITGB3",]>=1),union(which(pbmc5@assays$RNA@data["ITGA2B",]>=1),which(pbmc5@assays$RNA@data["PPBP",]>=1)))
pbmc5@meta.data$seurat_clusters[l]<-rep("platelet",length(l))

l<-union(which(pbmc5@assays$RNA@data["IGKC",]>=1),union(which(pbmc5@assays$RNA@data["IGHM",]>=1),which(pbmc5@assays$RNA@data["CD19",]>=1)))
pbmc5@meta.data$seurat_clusters[l]<-rep("Bcell",length(l))

l<-which(pbmc5@assays$RNA@data["CD14",]>=1)
pbmc5@meta.data$seurat_clusters[l]<-rep("CD14+mono",length(l))

l<-intersect(which(pbmc5@assays$RNA@data["FCGR3A",]>=0.5),which(pbmc5@meta.data$seurat_clusters==3))
pbmc5@meta.data$seurat_clusters[l]<-rep("CD16+mono",length(l))


l<-intersect(which(pbmc5@assays$RNA@data["CD3G",]<=0.3),intersect(which(pbmc5@assays$RNA@data["KLRF1",]>=0.5),which(pbmc5@meta.data$seurat_clusters %in% c("0","1","2"))))
pbmc5@meta.data$seurat_clusters[l]<-rep("NK",length(l))

l<-intersect(which(pbmc5@assays$RNA@data["CD8A",]<=0.3),intersect(which(pbmc5@assays$RNA@data["CD3G",]>0.5),intersect(which(pbmc5@assays$RNA@data["CD4",]>=0.5),which(pbmc5@meta.data$seurat_clusters %in% c("0","1","2")))))
pbmc5@meta.data$seurat_clusters[l]<-rep("CD4+T",length(l))

l<-intersect(which(pbmc5@assays$RNA@data["CD4",]<=0.3),intersect(which(pbmc5@assays$RNA@data["CD3G",]>0.5),intersect(which(pbmc5@assays$RNA@data["CD8A",]>=0.5),which(pbmc5@meta.data$seurat_clusters %in% c("0","1","2")))))
pbmc5@meta.data$seurat_clusters[l]<-rep("CD8+T",length(l))

l<-which(pbmc5$seurat_clusters %in% c("0","1","2","3"))
pbmc5@meta.data$seurat_clusters[l]<-rep("Unclear",length(l))
write.csv(pbmc5@meta.data,"PBMCgroup.csv")

FeaturePlot(pbmc5,"PPBP")
dim(pbmc5@meta.data)

pbmc6<-pbmc5
setwd("~\\PBMC")
group<-read.csv("PBMCgroup.csv")
identical(group[,1],rownames(pbmc6@meta.data))
pbmc6$celltype2<-pbmc6$seurat_clusters
pbmc6@active.ident<-as.factor(pbmc6$celltype2)
names(pbmc6@active.ident)<-rownames(pbmc6@meta.data)
UMAPPlot(pbmc6,pt.size=1.5)

#markers：
marker2gene<-c("CD19","CD14","FCGR3A","CD4","CD8A","KLRF1","PPBP")
marker2cell<-c("Bcell","CD14+mono","CD16+mono","CD4+T","CD8+T","NK","platelet")

library(reshape2)
vln.df=as.data.frame(pbmc6[["RNA"]]@data[marker2gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
pbmc6@meta.data$CB<-rownames(pbmc6@meta.data)
anno=pbmc6@meta.data[,c("CB","celltype2")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = marker2gene) #为了控制画图的基因顺序
vln.df$celltype=factor(vln.df$celltype,levels = marker2cell) #为了控制画图的基因顺序
ggplot(vln.df,aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(gene~.,scales = "free_y")+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 270,hjust = 1,vjust = 1),strip.text=element_text(size=35),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none",axis.text.x=element_text(size=50),axis.text.y=element_text(size=25),axis.title.y=element_text(size=50)
  )
ggsave("violinplot.png",width = 40,height = 65,limitsize=FALSE)

marker2gene<-c("CD19","CD14","FCGR3A","CD4","CD8A","KLRF1","PPBP","Unclear")
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
ggsave("cellnumber.png",gg,width = 22,height = 12,limitsize=FALSE)
gg<-ggplot(numberc,aes(x=celltype,y=ratio))+geom_bar(aes(fill=celltype),stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle = 45,size=25,vjust=0.5),axis.title.x=element_text(size=35),axis.text.y=element_text(size=35),axis.title.y=element_text(size=25),legend.position = "none")
ggsave("cellratio.png",gg,width = 22,height = 12,limitsize=FALSE)
