library(stringr)
library(ggplot2)
pastevector<-function(c,mark){
  if (length(c)==1)
  {
    k=c
  }
  if (length(c)>1)
  {
    k=c[1]
    for (i in c(2:length(c)))
    {
      k=paste(k,c[i],sep=mark)
    }
  }
  return(k)
} 

Generatecodon<-function(group=TRUE) 
{
  codons<-c()
  t=0
  for (i in c("T","C","A","G"))
  {
    for (j in c("T","C","A","G"))
    {
      for (k in c("T","C","A","G"))
      {
        if (pastevector(c(i,j,k),"") %in% c("TGA","TAG","TAA","ATG","TGG")==FALSE)
        {
          t=t+1
          codons[t]=pastevector(c(i,j,k),"")
        }
      }
    }
  }
  
  if (group==TRUE)
  {
    aas<-c("Phe","Phe","Leu1","Leu1","Ser2","Ser2","Ser2","Ser2","Tyr","Tyr","Cys","Cys","Leu2","Leu2","Leu2","Leu2","Pro","Pro","Pro","Pro","His","His","Gln","Gln","Arg2","Arg2","Arg2","Arg2","Ile","Ile","Ile","Thr","Thr","Thr","Thr","Asn","Asn","Lys","Lys","Ser1","Ser1","Arg1","Arg1","Val","Val","Val","Val","Ala","Ala","Ala","Ala","Asp","Asp","Glu","Glu","Gly","Gly","Gly","Gly")
  }
  if (group==FALSE)
  {
    aas<-c("Phe","Phe","Leu","Leu","Ser","Ser","Ser","Ser","Tyr","Tyr","Cys","Cys","Leu","Leu","Leu","Leu","Pro","Pro","Pro","Pro","His","His","Gln","Gln","Arg","Arg","Arg","Arg","Ile","Ile","Ile","Thr","Thr","Thr","Thr","Asn","Asn","Lys","Lys","Ser","Ser","Arg","Arg","Val","Val","Val","Val","Ala","Ala","Ala","Ala","Asp","Asp","Glu","Glu","Gly","Gly","Gly","Gly")
  }
  codonaas<-data.frame(aas)
  rownames(codonaas)<-codons
  colnames(codonaas)[1]<-"aas"
  return(codonaas)
}

RSCU<-function(highexpression,group=TRUE) 
{
  library(stringr)
  codontable<-Generatecodon(group)
  codontable$counts<-rep(0,dim(codontable)[1])
  codontable$RSCU<-rep(0,dim(codontable)[1])
  codontable$w<-rep(0,dim(codontable)[1])
  codontable$wadj<-rep(0,dim(codontable)[1])
  for (i in c(1:length(highexpression)))
  {
    for (j in c(1:(str_length(highexpression[i])/3)))
    {
      if (str_sub(highexpression[i],3*j-2,3*j) %in% c("TGA","TAG","TAA","ATG","TGG")==FALSE)
      {
        codontable[str_sub(highexpression[i],3*j-2,3*j),2]=codontable[str_sub(highexpression[i],3*j-2,3*j),2]+1
      }
    }
  }
  for (i in c(1:length(codontable[,1])))
  {
    codontable$RSCU[i]<-codontable$counts[i]/(mean(codontable$counts[which(codontable$aas==codontable$aas[i])]))
  }
  for (i in c(1:length(codontable[,1])))
  {
    codontable$w[i]<-codontable$RSCU[i]/max(codontable$RSCU[which(codontable$aas==codontable$aas[i])])
  }
  for  (i in c(1:length(codontable[,1])))
  {
    codontable$wadj[i]<-(codontable$counts[i]+0.1)/(max(codontable$counts[which(codontable$aas==codontable$aas[i])])+0.1*length(which(codontable$aas==codontable$aas[i])))
  }
  return(codontable)
}


CAI<-function(sequence,codontable,adj=TRUE) 
{
  library(stringr)
  sequence=str_to_upper(str_replace_all(str_to_upper(sequence),"U","T"))
  n=4
  if (adj==TRUE)
  {
    n=5
  }
  logCAI=0
  l=0
  for (j in c(1:(str_length(sequence)/3)))
  {
    if ((str_sub(sequence,3*j-2,3*j) %in% c("TGA","TAG","TAA","ATG","TGG")==FALSE) && (str_sub(sequence,3*j-2,3*j) %in% rownames(codontable)))
    {
      logCAI=logCAI+log10(codontable[str_sub(sequence,3*j-2,3*j),n])
      l=l+1
    }
  }
  CAI=10^(logCAI/l)
  return(CAI)
}

###SARS-COV-2
samples<-c("Control4h","Control96h","SARS4h","SARS24h","SARS48h","SARS72h","SARS96h")
for (i in c(1:length(samples)))
{
  setwd(paste("~/RiboseqSARS/RNA",samples[i],sep="/"))
  TOP200<-read.csv(paste(samples[i],"top200codinggenemaxtranscriptCDS.csv",sep=""))
  TOP5000<-read.csv(paste(samples[i],"MajorCDSforalltranscripts(FPKM>1).csv",sep=""))
  codontable<-RSCU(TOP200[,5],group=TRUE)
  TOP5000CAI<-data.frame(gene=TOP5000$gene,seq=TOP5000$CDSsequence,CAI=rep(0,length(TOP5000$gene)))
  for (j in c(1:dim(TOP5000CAI)[1]))
  {
    TOP5000CAI$CAI[j]<-CAI(TOP5000CAI$seq[j],codontable,adj=TRUE)
  }
  write.csv(TOP5000CAI,paste(samples[i],"TOP5000genetranscriptCAI.csv",sep=""))
}

###HIV
samples<-c("Control","HIV")
for (i in c(1:length(samples)))
{
  setwd(paste("~/RiboseqHIV/RNA",samples[i],sep="/"))
  TOP200<-read.csv(paste(samples[i],"top200codinggenemaxtranscriptCDS.csv",sep=""))
  TOP5000<-read.csv(paste(samples[i],"MajorCDSforalltranscripts(FPKM>1).csv",sep=""))
  codontable<-RSCU(TOP200[,5],group=TRUE)
  TOP5000CAI<-data.frame(gene=TOP5000$gene,seq=TOP5000$CDSsequence,CAI=rep(0,length(TOP5000$gene)))
  for (j in c(1:dim(TOP5000CAI)[1]))
  {
    TOP5000CAI$CAI[j]<-CAI(TOP5000CAI$seq[j],codontable,adj=TRUE)
  }
  write.csv(TOP5000CAI,paste(samples[i],"TOP5000genetranscriptCAI.csv",sep=""))
}

library(rtracklayer)
setwd("~")
gtf<-import("gencode.v39.annotation.gtf")
gtf<-as.data.frame(gtf)
#SARS-CoV-2
samples<-c("Control4h","Control96h","SARS4h","SARS24h","SARS48h","SARS72h","SARS96h")
for (i in c(1:length(samples)))
{
  setwd(paste("~/RiboseqSARS/RNA",samples[i],sep="/"))
  RiboFPKM<-read.csv("RiboFPKM.csv")
  RNAFPKM<-read.csv("RNAFPKM.csv")
  TOP5000CAI<-read.csv(paste(samples[i],"TOP5000genetranscriptCAI.csv",sep=""))
  TOP5000CAI$TE<-rep(0,dim(TOP5000CAI)[1])
  TOP5000CAI$gene<-gtf$gene_name[match(TOP5000CAI$gene,gtf$gene_id)]
  TOP5000CAI$RNAFPKM<-RNAFPKM$meanFPKM[match(TOP5000CAI$gene,RNAFPKM$geneid)]
  TOP5000CAI$RiboFPKM<-RiboFPKM$meanFPKM[match(TOP5000CAI$gene,RiboFPKM$geneid)]
  TOP5000CAI$TE<-TOP5000CAI$RiboFPKM/TOP5000CAI$RNAFPKM
  TOP5000CAI$logTE<-log(TOP5000CAI$TE)
  #write.csv(TOP5000CAI,paste(samples[i],"TOP5000CAIandTE.csv",sep=""))
  TECAImodel<-lm(logTE~CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),])
  g<-ggplot(TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],aes(x=CAI,y=logTE))+geom_point()+theme_bw()+theme(axis.title = element_text(size=20),axis.text=element_text(size=18))+ geom_smooth(method='lm', formula= y~x)+annotate("text",x=0.9,y=5,size=5,label=paste(paste(paste(paste("logTE=","*CAI",sep=as.character(round(unname(TECAImodel$coefficients[2]),3))),as.character(abs(round(unname(TECAImodel$coefficients[1]),3))),sep=(if (TECAImodel$coefficients[1]>0) {"+"} else {"-"})),paste("rho=",round(cor.test(~logTE+CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],method="spearman")$estimate,3),sep=""),sep="\n"),paste("p=",signif(cor.test(~logTE+CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],method="spearman")$p.value,3),sep=""),sep=","),colour="red")+scale_x_continuous(limits =c(0.5,1.05),breaks=c(0.5,0.6,0.7,0.8,0.9,1.0),labels=c("0.5","0.6","0.7","0.8","0.9","1.0"),expand=c(0,0))
  ggsave(paste(samples[i],"CAIandTEregression.png",sep=""),g,width=8,height=8)
}
#HIV-1
samples<-c("Control","HIV")
for (i in c(1:length(samples)))
{
  setwd(paste("~/RiboseqHIV/RNA",samples[i],sep="/"))
  RiboFPKM<-read.csv("RiboFPKM.csv")
  RNAFPKM<-read.csv("RNAFPKM.csv")
  TOP5000CAI<-read.csv(paste(samples[i],"TOP5000genetranscriptCAI.csv",sep=""))
  TOP5000CAI$TE<-rep(0,dim(TOP5000CAI)[1])
  TOP5000CAI$gene<-gtf$gene_name[match(TOP5000CAI$gene,gtf$gene_id)]
  TOP5000CAI$RNAFPKM<-RNAFPKM$meanFPKM[match(TOP5000CAI$gene,RNAFPKM$geneid)]
  TOP5000CAI$RiboFPKM<-RiboFPKM$meanFPKM[match(TOP5000CAI$gene,RiboFPKM$geneid)]
  TOP5000CAI$TE<-TOP5000CAI$RiboFPKM/TOP5000CAI$RNAFPKM
  TOP5000CAI$logTE<-log(TOP5000CAI$TE)
  #write.csv(TOP5000CAI,paste(samples[i],"TOP5000CAIandTE.csv",sep=""))
  TECAImodel<-lm(logTE~CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),])
  g<-ggplot(TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],aes(x=CAI,y=logTE))+geom_point()+theme_bw()+theme(axis.title = element_text(size=20),axis.text=element_text(size=18))+ geom_smooth(method='lm', formula= y~x)+annotate("text",x=0.9,y=5,size=5,label=paste(paste(paste(paste("logTE=","*CAI",sep=as.character(round(unname(TECAImodel$coefficients[2]),3))),as.character(abs(round(unname(TECAImodel$coefficients[1]),3))),sep=(if (TECAImodel$coefficients[1]>0) {"+"} else {"-"})),paste("rho=",round(cor.test(~logTE+CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],method="spearman")$estimate,3),sep=""),sep="\n"),paste("p=",signif(cor.test(~logTE+CAI,data=TOP5000CAI[which(is.finite(TOP5000CAI$logTE)==TRUE),],method="spearman")$p.value,3),sep=""),sep=","),colour="red")+scale_x_continuous(limits =c(0.5,1.05),breaks=c(0.5,0.6,0.7,0.8,0.9,1.0),labels=c("0.5","0.6","0.7","0.8","0.9","1.0"),expand=c(0,0))
  ggsave(paste(samples[i],"CAIandTEregression.png",sep=""),g,width=8,height=8)
}