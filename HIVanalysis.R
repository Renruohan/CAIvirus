library(ggplot2)
library(clusterProfiler)
library(stringr)
setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()
  
gag<-c()
gagpol<-c()
env<-c()
vif<-c()
vpr<-c()
tat<-c()
rev<-c()
vpu<-c()
nef<-c()
for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    gag<-c(gag,cviruscai$CAIadj[1])
    gagpol<-c(gagpol,cviruscai$CAIadj[2])
    env<-c(env,cviruscai$CAIadj[3])
    vif<-c(vif,cviruscai$CAIadj[4])
    vpr<-c(vpr,cviruscai$CAIadj[5])
    tat<-c(tat,cviruscai$CAIadj[6])
    rev<-c(rev,cviruscai$CAIadj[7])
    vpu<-c(vpu,cviruscai$CAIadj[8])
    nef<-c(nef,cviruscai$CAIadj[9])
  }
}

gag<-gag[intersect(which(is.na(gag)==FALSE),which(gag!="inapplicability"))]
gagpol<-gagpol[intersect(which(is.na(gagpol)==FALSE),which(gagpol!="inapplicability"))]
env<-env[intersect(which(is.na(env)==FALSE),which(env!="inapplicability"))]
vif<-vif[intersect(which(is.na(vif)==FALSE),which(vif!="inapplicability"))]
vpr<-vpr[intersect(which(is.na(vpr)==FALSE),which(vpr!="inapplicability"))]
tat<-tat[intersect(which(is.na(tat)==FALSE),which(tat!="inapplicability"))]
rev<-rev[intersect(which(is.na(rev)==FALSE),which(rev!="inapplicability"))]
vpu<-vpu[intersect(which(is.na(vpu)==FALSE),which(vpu!="inapplicability"))]
nef<-nef[intersect(which(is.na(nef)==FALSE),which(nef!="inapplicability"))]

gene<-c(rep("gag",length(gag)),rep("gagpol",length(gagpol)),rep("env",length(env)),rep("vif",length(vif)),rep("vpr",length(vpr)),
rep("tat",length(tat)),rep("rev",length(rev)),rep("vpu",length(vpu)),rep("nef",length(nef)))
cai<-c(gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)
caidata<-data.frame(gene,as.numeric(cai))
colnames(caidata)<-c("Gene","CAI")
gg<-ggplot(caidata,aes(x=Gene,y=CAI))+geom_violin(aes(fill=Gene))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=28),legend.title=element_text(size=30))
ggsave("HIVallgeneallCAI.png",gg, width = 20, height = 8)


mean(as.numeric(gag))
mean(as.numeric(gagpol))
mean(as.numeric(env))
mean(as.numeric(vif))
mean(as.numeric(vpr))
mean(as.numeric(tat))
mean(as.numeric(rev))
mean(as.numeric(vpu))
mean(as.numeric(nef))
min(caidata[,2])
max(caidata[,2])


setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()

for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  gag<-c()
  gagpol<-c()
  env<-c()
  vif<-c()
  vpr<-c()
  tat<-c()
  rev<-c()
  vpu<-c()
  nef<-c()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    gag<-c(gag,cviruscai$CAIadj[1])
    gagpol<-c(gagpol,cviruscai$CAIadj[2])
    env<-c(env,cviruscai$CAIadj[3])
    vif<-c(vif,cviruscai$CAIadj[4])
    vpr<-c(vpr,cviruscai$CAIadj[5])
    tat<-c(tat,cviruscai$CAIadj[6])
    rev<-c(rev,cviruscai$CAIadj[7])
    vpu<-c(vpu,cviruscai$CAIadj[8])
    nef<-c(nef,cviruscai$CAIadj[9])
  }
  gag<-gag[intersect(which(is.na(gag)==FALSE),which(gag!="inapplicability"))]
  gagpol<-gagpol[intersect(which(is.na(gagpol)==FALSE),which(gagpol!="inapplicability"))]
  env<-env[intersect(which(is.na(env)==FALSE),which(env!="inapplicability"))]
  vif<-vif[intersect(which(is.na(vif)==FALSE),which(vif!="inapplicability"))]
  vpr<-vpr[intersect(which(is.na(vpr)==FALSE),which(vpr!="inapplicability"))]
  tat<-tat[intersect(which(is.na(tat)==FALSE),which(tat!="inapplicability"))]
  rev<-rev[intersect(which(is.na(rev)==FALSE),which(rev!="inapplicability"))]
  vpu<-vpu[intersect(which(is.na(vpu)==FALSE),which(vpu!="inapplicability"))]
  nef<-nef[intersect(which(is.na(nef)==FALSE),which(nef!="inapplicability"))]
  gene<-c(rep("gag",length(gag)),rep("gagpol",length(gagpol)),rep("env",length(env)),rep("vif",length(vif)),rep("vpr",length(vpr)),
          rep("tat",length(tat)),rep("rev",length(rev)),rep("vpu",length(vpu)),rep("nef",length(nef)))
  cai<-c(gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)
  caidata<-data.frame(gene,as.numeric(cai))
  colnames(caidata)<-c("Gene","CAI")
  gg<-ggplot(caidata,aes(x=Gene,y=CAI))+geom_violin(aes(fill=Gene))+theme_classic()+theme(axis.text=element_text(colour='black',size=20),axis.title=element_text(colour='black',size=20))
  setwd("~\\CAI\\HIVanalysis")
  gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=28),legend.title=element_text(size=30))
  ggsave(paste("HIVallgene","CAI.png",sep=cell),gg, width = 20, height = 6)
}




library(ggplot2)
library(clusterProfiler)
setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()

gag<-c()
gagpol<-c()
env<-c()
vif<-c()
vpr<-c()
tat<-c()
rev<-c()
vpu<-c()
nef<-c()
cells<-c()
for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    gag<-c(gag,cviruscai$CAIadj[1])
    gagpol<-c(gagpol,cviruscai$CAIadj[2])
    env<-c(env,cviruscai$CAIadj[3])
    vif<-c(vif,cviruscai$CAIadj[4])
    vpr<-c(vpr,cviruscai$CAIadj[5])
    tat<-c(tat,cviruscai$CAIadj[6])
    rev<-c(rev,cviruscai$CAIadj[7])
    vpu<-c(vpu,cviruscai$CAIadj[8])
    nef<-c(nef,cviruscai$CAIadj[9])
    cells<-c(cells,cell)
  }
}

gag<-data.frame(gag,cells)
gagpol<-data.frame(gagpol,cells)
env<-data.frame(env,cells)
vif<-data.frame(vif,cells)
vpr<-data.frame(vpr,cells)
tat<-data.frame(tat,cells)
rev<-data.frame(rev,cells)
vpu<-data.frame(vpu,cells)
nef<-data.frame(nef,cells)



gag<-gag[intersect(which(is.na(gag[,1])==FALSE),which(gag[,1]!="inapplicability")),]
gag[,3]<-rep("gag",dim(gag)[1])
gagpol<-gagpol[intersect(which(is.na(gagpol[,1])==FALSE),which(gagpol[,1]!="inapplicability")),]
gagpol[,3]<-rep("gagpol",dim(gagpol)[1])
env<-env[intersect(which(is.na(env[,1])==FALSE),which(env[,1]!="inapplicability")),]
env[,3]<-rep("env",dim(env)[1])
vif<-vif[intersect(which(is.na(vif[,1])==FALSE),which(vif[,1]!="inapplicability")),]
vif[,3]<-rep("vif",dim(vif)[1])
vpr<-vpr[intersect(which(is.na(vpr[,1])==FALSE),which(vpr[,1]!="inapplicability")),]
vpr[,3]<-rep("vpr",dim(vpr)[1])
tat<-tat[intersect(which(is.na(tat[,1])==FALSE),which(tat[,1]!="inapplicability")),]
tat[,3]<-rep("tat",dim(tat)[1])
rev<-rev[intersect(which(is.na(rev[,1])==FALSE),which(rev[,1]!="inapplicability")),]
rev[,3]<-rep("rev",dim(rev)[1])
vpu<-vpu[intersect(which(is.na(vpu[,1])==FALSE),which(vpu[,1]!="inapplicability")),]
vpu[,3]<-rep("vpu",dim(vpu)[1])
nef<-nef[intersect(which(is.na(nef[,1])==FALSE),which(nef[,1]!="inapplicability")),]
nef[,3]<-rep("nef",dim(nef)[1])
colnames(gag)<-c("CAI","Celltype","Gene")
colnames(gagpol)<-c("CAI","Celltype","Gene")
colnames(env)<-c("CAI","Celltype","Gene")
colnames(vif)<-c("CAI","Celltype","Gene")
colnames(vpr)<-c("CAI","Celltype","Gene")
colnames(tat)<-c("CAI","Celltype","Gene")
colnames(rev)<-c("CAI","Celltype","Gene")
colnames(vpu)<-c("CAI","Celltype","Gene")
colnames(nef)<-c("CAI","Celltype","Gene")
c<-rbind(gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)

c1<-c[grep(c$Celltype,pattern="lung"),]
gg<-ggplot(c1,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgenelunggroupingCAI.png",gg, width = 30, height = 9,limitsize=FALSE)

c1<-c[grep(c$Celltype,pattern="PBMC"),]
gg<-ggplot(c1,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgenePBMCgroupingCAI.png",gg, width = 30, height = 8,limitsize=FALSE)

c1<-c[grep(c$Celltype,pattern="CD4"),]
c1<-c1[which(str_detect(c1$Celltype,"CD4\\+")==FALSE),]
gg<-ggplot(c1,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgeneCD4groupingCAI.png",gg, width = 14, height = 8,limitsize=FALSE)

c1<-c[which(c$Celltype %in% c("bloodmonocyte","dermalmacro","intestinemacro","kupffercell","l-dermalDC","l+dermalDC","langerhanscell" ,"microgoliaCC","microgoliaOC","monoderivedDC", "monoderivedmacro","monoderivedosteoclast")),]
gg<-ggplot(c1,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgenemyeloidgroupingCAI.png",gg, width = 20, height = 8,limitsize=FALSE)


gg<-ggplot(c,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgenegroupingCAI.png",gg, width = 100, height = 8,limitsize=FALSE)

c1<-c[which(c$Celltype %in% c("adipocyte","hepatocyte","cholangiocyteCD133+","cholangiocyteCD133-","satelite")),]
c1$Celltype[which(c1$Celltype=="cholangiocyteCD133+")]<-rep("CD133+ cholangiocyte",length(c1$Celltype[which(c1$Celltype=="cholangiocyteCD133+")]))
c1$Celltype[which(c1$Celltype=="cholangiocyteCD133-")]<-rep("CD133- cholangiocyte",length(c1$Celltype[which(c1$Celltype=="cholangiocyteCD133-")]))
c1$Celltype[which(c1$Celltype=="satelite")]<-rep("hepatic satellite cell",length(c1$Celltype[which(c1$Celltype=="satelite")]))
gg<-ggplot(c1,aes(x=Gene,y=as.numeric(CAI)))+geom_violin(aes(fill=Celltype))+theme_classic()+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18))
setwd("~\\CAI\\HIVanalysis")
gg<-gg+xlab("ORF")+ylab("CAI")+theme(axis.text=element_text(colour='black',size=18),axis.title=element_text(colour='black',size=18),axis.text.x=element_text(size=28),axis.text.y=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28),legend.text=element_text(size=24),legend.title=element_text(size=30))
ggsave("HIVallgenemetabolismgroupingCAI.png",gg, width = 14, height = 8,limitsize=FALSE)



setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()
gag<-c()
gagpol<-c()
env<-c()
vif<-c()
vpr<-c()
tat<-c()
rev<-c()
vpu<-c()
nef<-c()
virusname<-c()
for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    if (("inapplicability" %in% cviruscai$CAIadj==FALSE) && (length(which(is.na(cviruscai$CAIadj)))==0))
    {
      gag<-c(gag,cviruscai$CAIadj[1])
      gagpol<-c(gagpol,cviruscai$CAIadj[2])
      env<-c(env,cviruscai$CAIadj[3])
      vif<-c(vif,cviruscai$CAIadj[4])
      vpr<-c(vpr,cviruscai$CAIadj[5])
      tat<-c(tat,cviruscai$CAIadj[6])
      rev<-c(rev,cviruscai$CAIadj[7])
      vpu<-c(vpu,cviruscai$CAIadj[8])
      nef<-c(nef,cviruscai$CAIadj[9])
      virusname<-c(virusname,virus[j])
    }
  }
}

allCDS<-data.frame(gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)
corallCDS <- cor(allCDS)

library(tidyverse)
library(reshape2)
data <- as.data.frame(corallCDS) %>% 
  mutate(x=rownames(corallCDS)) %>% 
  melt(id='x') %>%                 
  rename('y'='variable','Corr'='value')


list <- rownames(corallCDS)
list <- factor(list,levels = list)

ggplot(data,aes(factor(x,levels = list),
                factor(y,levels = list),
                fill=Corr))+  
  geom_tile()+
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red',
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL)+
  theme_bw(base_size = 15)





#evolution
setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()
gag<-c()
gagpol<-c()
env<-c()
vif<-c()
vpr<-c()
tat<-c()
rev<-c()
vpu<-c()
nef<-c()
allCDS<-data.frame()
virusname<-c()
for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  gag<-c()
  gagpol<-c()
  env<-c()
  vif<-c()
  vpr<-c()
  tat<-c()
  rev<-c()
  vpu<-c()
  nef<-c()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    if (("inapplicability" %in% cviruscai$CAIadj==FALSE) && (length(which(is.na(cviruscai$CAIadj)))==0))
    {
    gag<-c(gag,cviruscai$CAIadj[1])
    gagpol<-c(gagpol,cviruscai$CAIadj[2])
    env<-c(env,cviruscai$CAIadj[3])
    vif<-c(vif,cviruscai$CAIadj[4])
    vpr<-c(vpr,cviruscai$CAIadj[5])
    tat<-c(tat,cviruscai$CAIadj[6])
    rev<-c(rev,cviruscai$CAIadj[7])
    vpu<-c(vpu,cviruscai$CAIadj[8])
    nef<-c(nef,cviruscai$CAIadj[9])
    if (i==1)
    {
       virusname<-c(virusname,virus[j])
       allCDS<-data.frame(virusname)
    }
    }
  }
  allCDS<-data.frame(allCDS,gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)
}

setwd("~\\CAI")
HIVvar<-read.csv("HIVvariantsequence.csv")
setwd("~\\CAI\\HIVsequence")

library(tidyverse)
library(reshape2)
rownames(allCDS)<-allCDS[,1]
cleandata<-allCDS[,-1]
data.pca <- prcomp(cleandata, center = F, scale. = F)
dd<-as.data.frame(data.pca$x)
dd$sample<-rownames(dd)
dd
Accessions<-str_sub(rownames(dd),1,8)
HIVsubtype<-data.frame(HIVvar$Accession[match(Accessions,HIVvar$Accession)],HIVvar$Subtype[match(Accessions,HIVvar$Accession)])
identical(str_sub(rownames(dd),1,8),HIVsubtype[,1])
dd$subtype<-HIVsubtype[,2]
ggplot(dd,aes(x=PC1,y=PC2)) + geom_point(aes(col=subtype),size=5)+theme_classic()
rownames(cleandata)<-str_sub(rownames(dd),1,8)


library(ggplot2)
library(clusterProfiler)
library(stringr)
library(ggheatmap)
setwd("~\\CAI\\HIVsequenceCDS")
celltypes<-list.files()
gag<-c()
gagpol<-c()
env<-c()
vif<-c()
vpr<-c()
tat<-c()
rev<-c()
vpu<-c()
nef<-c()
allCDS<-data.frame()
virusname<-c()
for (i in c(1:length(celltypes)))
{
  cell<-celltypes[i]
  setwd(paste("~\\CAI\\HIVsequenceCDS",cell,sep="\\"))
  virus<-list.files()
  gag<-c()
  gagpol<-c()
  env<-c()
  vif<-c()
  vpr<-c()
  tat<-c()
  rev<-c()
  vpu<-c()
  nef<-c()
  for (j in c(1:length(virus)))
  {
    cviruscai<-read.csv(virus[j])
    if (("inapplicability" %in% cviruscai$CAIadj==FALSE) && (length(which(is.na(cviruscai$CAIadj)))==0))
    {
      gag<-c(gag,cviruscai$CAIadj[1])
      gagpol<-c(gagpol,cviruscai$CAIadj[2])
      env<-c(env,cviruscai$CAIadj[3])
      vif<-c(vif,cviruscai$CAIadj[4])
      vpr<-c(vpr,cviruscai$CAIadj[5])
      tat<-c(tat,cviruscai$CAIadj[6])
      rev<-c(rev,cviruscai$CAIadj[7])
      vpu<-c(vpu,cviruscai$CAIadj[8])
      nef<-c(nef,cviruscai$CAIadj[9])
      if (i==1)
      {
        virusname<-c(virusname,virus[j])
        allCDS<-data.frame(virusname)
      }
    }
  }
  allCDS<-data.frame(allCDS,gag,gagpol,env,vif,vpr,tat,rev,vpu,nef)
}
rownames(allCDS)<-allCDS[,1]
Accessions<-str_sub(rownames(allCDS),1,8)
setwd("~\\CAI")
HIVvar<-read.csv("HIVvariantsequence.csv")
setwd("~\\CAI\\HIVsequence")
HIVsubtype<-data.frame(HIVvar$Accession[match(Accessions,HIVvar$Accession)],HIVvar$Subtype[match(Accessions,HIVvar$Accession)])
identical(str_sub(rownames(allCDS),1,8),HIVsubtype[,1])
rownames(allCDS)<-HIVsubtype[,1]
library(tidyverse)
library(reshape2)

cleandata<-allCDS[,-1]
cleandata<-t(cleandata)
dd10<-cleandata
annotation_col = data.frame(HIVsubtype[,2])
rownames(annotation_col) = colnames(dd10)
colnames(annotation_col)<-"Strain"
annotation_row<-data.frame(genes=rownames(dd10))
rownames(annotation_row)<-rownames(dd10)


c<- list(Strain=c("C" = "#0000FF","A1" = "#99CC32","B"="#5C3317","L"="#D9D919","N"="#70DBDB","G"="#00FF7F","O"="#DB9370","D"="#D8D8BF","H"="#70DB93","F1"="#FF0000","J"="#FF00FF","A6"="#38B0DE","F2"="#00009C","A7"="#CDCDCD","U"="#EBC79E","A2"="#FFFFFF"))
        
setwd("~\\CAI\\HIVanalysis")
dd11<-apply(dd10,2,function(x){(x-mean(x))/sd(x)})
p<-ggheatmap(dd11,cluster_cols = T,levels_cols =colnames(dd10),levels_rows = rev(rownames(dd10)),color = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),annotation_cols = annotation_col,annotation_color = c) %>% ggheatmap_theme(1,theme=list(theme(axis.text.x = element_text(angle = 72,hjust = 1,vjust = 1,color="black", size=12),legend.title=element_text(size=20),legend.text=element_text(size=15),axis.text.y = element_text(color="black", size=0.1))))
p<-p %>% ggheatmap_theme(2,theme=list(theme(legend.title=element_text(size=20),legend.text=element_text(size=15),axis.text.y = element_text(color="black", size=0.1))))
ggsave("evolutionbyCAI.tiff",p,limitsize=FALSE,dpi=300,width=14,height=20)

