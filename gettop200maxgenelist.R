#scRNA-seq datasets of PBMC as example

library(stringr)
library(dplyr)
setwd("//home//haoyu//Desktop//PBMC")
class<-read.csv("PBMCgroup.csv")
celltype<-class$seurat_clusters
celltype<-celltype[!duplicated(celltype)]
setwd("//home//haoyu//Desktop")
gtf_data<-read.csv("gtf_dataframe.csv",row.names=1) 
gtfcodingtrans<-gtf_data[which(gtf_data$transcript_type=="protein_coding"),]
gtfCDS<-gtf_data[which(gtf_data$type=="CDS"),]
compare<-read.csv("gtfcompare.csv")
comparegene<-compare[which(compare$type=="gene"),] 
comparecodinggene<-comparegene[which(comparegene$gene_type=="protein_coding"),] 
setwd("//home//haoyu//Desktop//PBMC")

for (i in c(1:length(celltype)))
{
  setwd("//home//haoyu//Desktop//PBMC")
  ijt=class[which(class$seurat_clusters==celltype[i]),]
  OC<-read.table(ijt$X[1],header=TRUE)[,1:6]
  for (j in c(1:dim(ijt)[1]))
  {
    OC[,6+j]<-read.table(ijt$X[j],header=TRUE)[,7]
  }
  colnames(OC)[1:6]<-c("geneid","chr","start","end","strand","length")
  OCFPKM<-OC
  t<-rep(0,dim(ijt)[1])
  for(j in c(1:dim(ijt)[1]))
  {
    t[j]<-sum(OC[,6+j])
  }
  OCFPKM$meanFPKM<-rep(0,dim(OCFPKM)[1])
  OCFPKM$varFPKM<-rep(0,dim(OCFPKM)[1])
  OCFPKM$ID<-rep("0",dim(OCFPKM)[1])
  OCFPKM$coding<-rep("no",dim(OCFPKM)[1])
  for (j in c(1:length(OC[,1])))
  {
    OCFPKM[j,7:(dim(OCFPKM)[2]-4)]<-as.numeric(OC[j,7:(dim(OCFPKM)[2]-4)])*1000000000/(t*as.numeric(OC$length[j]))
    OCFPKM$meanFPKM[j]<-sum(OCFPKM[j,7:(dim(OCFPKM)[2]-4)])/dim(OCFPKM[j,7:(dim(OCFPKM)[2]-4)])[2]
    OCFPKM$varFPKM[j]<-sd(OCFPKM[j,7:(dim(OCFPKM)[2]-4)])
    if (OCFPKM$geneid[j] %in% comparegene$gene_name)
    {
      OCFPKM$ID[j]<-comparegene$gene_id[which(comparegene$gene_name==OCFPKM$geneid[j])][1]
    }
    if (OCFPKM$geneid[j] %in% comparecodinggene$gene_name)
    {
      OCFPKM$coding[j]<-"yes"
    }
  }
  write.csv(OC,paste("human","counts.csv",sep=celltype[i]))
  write.csv(OCFPKM,paste("human","FPKM.csv",sep=celltype[i]))
  length(OCFPKM$ID[which(OCFPKM$ID!="0")])
  length(OCFPKM$ID[which(OCFPKM$coding!="no")])
  ###OCFPKM1<-OCFPKM[which(OCFPKM$varFPKM < 0.5*OCFPKM$meanFPKM),] #remove genes with too high variation
  OCFPKM1<-OCFPKM
  OCFPKM1<-OCFPKM1[which(OCFPKM1$coding=="yes"),] #remove non-coding genes
  OCFPKM1<-OCFPKM1[which(OCFPKM1$ID %in% gtfcodingtrans$gene_id),] 
  OCFPKM1<-OCFPKM1[which(OCFPKM1$ID %in% gtfCDS$gene_id),]
  OCFPKM1$ccdsid<-rep(0,dim(OCFPKM1)[1])
  for (j in c(1:dim(OCFPKM1)[1]))
  {
    tt=gtfCDS$ccdsid[which(gtfCDS$gene_id==OCFPKM1$ID[j])]
    if (length(which(is.na(tt)))<length(tt))
    {
      OCFPKM1$ccdsid[j]<-1 ##remove low-quality genes without annotation
    }
  }
  OCFPKM1<-OCFPKM1[which(OCFPKM1$ccdsid==1),]
  dim(OCFPKM1)
  maxgenefpkm<-arrange(OCFPKM1,desc(meanFPKM))[1:200,]
  setwd("//home//haoyu//Desktop//PBMC//maxgene")
  write.csv(maxgenefpkm,paste("human","maxgene",sep=celltype[i]))
}












