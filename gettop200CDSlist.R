#scRNA-seq datasets of PBMC as example

library(stringr)
library(dplyr)
library(stringr)
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
reverse<-function(sequence,strand)
{
  table<-data.frame(c("A","G","C","T","N"))
  rownames(table)<-c("T","C","G","A","N")
  library(stringr)
  if (strand=="+")
  {
    c=sequence
  }
  if (strand=="-")
  {
    l=str_length(sequence)
    c=sequence
    for (i in c(l:1))
    {
      str_sub(c,l+1-i,l+1-i)=table[str_sub(sequence,i,i),]
    }
  }
  return(c)
} 
setwd("//home//haoyu//Desktop//PBMC")
class<-read.csv("PBMCgroup.csv")
celltype<-class$seurat_clusters
celltype<-celltype[!duplicated(celltype)]
setwd("//home//haoyu//Desktop")
compare<-read.csv("gtfsplicingcompare.csv")
compare<-compare[which(compare$type=="exon"),] 
gene_idd<-read.csv("gene_idd.csv",row.names=1)[,1]
tosingle<-read.csv("tosingle.csv",row.names=1)[,1]
compare<-compare[tosingle,]
compare$gene_id<-gene_idd
gtf_data<-read.csv("gtf_dataframe.csv",row.names=1) 
gtfcodingtrans<-gtf_data[which(gtf_data$transcript_type=="protein_coding"),]
gtfCDS<-gtf_data[which(gtf_data$type=="CDS"),]
chromosome<-data.frame(c(1:25))
rownames(chromosome)<-as.character(gtfCDS$seqnames[!duplicated(gtfCDS$seqnames)])
a<-read.csv("gtfcompare.csv")
fasta<-readr::read_csv("GRCh38.p13.genome.fa.csv") 



for (i in c(1:length(celltype)))
{
  setwd("//home//haoyu//Desktop//PBMC//maxgene")
  top200<-read.csv(paste("human","maxgene",sep=celltype[i]))
  ijt=class[which(class$seurat_clusters==celltype[i]),]
  paste(str_split(ijt$X[1],"\\.")[[1]][1],"out",sep=".AS.")
  setwd("//home//haoyu//Desktop//PBMC//AS")
  OC<-read.table(paste(str_split(ijt$X[1],"\\.")[[1]][1],"out",sep=".AS."),header=TRUE)[,1:6]
  for (j in c(1:dim(ijt)[1]))
  {
    OC[,6+j]<-read.table(paste(str_split(ijt$X[j],"\\.")[[1]][1],"out",sep=".AS."),header=TRUE)[,7]
  }
  colnames(OC)[1:6]<-c("geneid","chr","start","end","strand","length")
  t<-rep(0,dim(ijt)[1])
  for(j in c(1:dim(ijt)[1]))
  {
    t[j]<-sum(OC[,6+j])
  }
  OC<-OC[tosingle,]
  identical(compare$gene_id[1:10000],OC$geneid[1:10000]) #FALSE
  OC$geneid<-gene_idd
  identical(compare$gene_id,OC$geneid)
  sum(str_count(compare$gene_id,"\\+")) #0
  
  
  OC<-OC[which(OC$geneid %in% top200$ID),] 
  length(OC$geneid[!duplicated(OC$geneid)]) 
  comparetop200exon<-compare[which(compare$gene_id %in% top200$ID),] #top200gene annotation
  length(comparetop200exon$gene_id[!duplicated(comparetop200exon$gene_id)]) #200
  transcriptstop200<-c() 
  for (j in c(1:length(comparetop200exon[,1])))
  {
    transcriptstop200<-c(transcriptstop200,str_split(comparetop200exon$transcripts[j],"\\+")[[1]])
  }
  transcriptstop200<-transcriptstop200[!duplicated(transcriptstop200)]
  finaltrans<-c()
  finalgene<-c()
  for (j in c(1:length(transcriptstop200)))
  {
    dup<-comparetop200exon$gene_id[which(str_detect(comparetop200exon$transcripts,transcriptstop200[j]))]
    dup<-dup[!duplicated(dup)]
    finaltrans<-c(finaltrans,rep(transcriptstop200[j],length(dup)))
    finalgene<-c(finalgene,dup)
  }
  finallist<-data.frame(finaltrans,finalgene)
  colnames(finallist)<-c("transcript_id","gene_id")
  identical(comparetop200exon$start,OC$start)
  identical(comparetop200exon$end,OC$end) 
  identical(comparetop200exon$strand,OC$strand) #TRUE
  length(finallist$gene_id[!duplicated(finallist$gene_id)]) #200
  
  
  
  
  
  
  
  
  
  ####OCFPKM[j,7:(dim(OC)[2]-1)]
  
  OC$transcripts<-comparetop200exon$transcripts
  colnames(OC)[c(7:(dim(OC)[2]-1))]<-paste("OC",c(1:(dim(OC)[2]-7)),sep="")
  
  for (j in c(7:(dim(OC)[2]-1)))
  {
    finallist[,j-4]<-rep(0,dim(finallist)[1])
    colnames(finallist)[j-4]<-paste("OC",j-6,sep="")
  }
  finallist$meancounts<-rep(0,dim(finallist)[1])
  for (j in c(7:(dim(OC)[2]-1))) 
  {
    finallist[,j-3+length(c(7:(dim(OC)[2]-1)))]<-rep(0,dim(finallist)[1])
    colnames(finallist)[j-3+length(c(7:(dim(OC)[2]-1)))]<-paste("OCFPKM",j-6,sep="")
  }
  finallist$meanFPKM<-rep(0,dim(finallist)[1])
  finallist$length<-rep(0,dim(finallist)[1])
  
  
  for (tt in c(1:length(finallist[,1])))
  {
    for (j in c(7:(dim(OC)[2]-1)))
    {
      finallist[tt,j-4]<-sum(OC[which(str_detect(OC$transcripts,finallist$transcript_id[tt])),j])
    }
    finallist$meancounts[tt]<-mean(as.numeric(finallist[tt,c(3:(dim(OC)[2]-5))]))
  }
  
  for (tt in c(1:length(finallist[,1])))
  {
    finallist$length[tt]<-sum(OC$length[which(str_detect(OC$transcripts,finallist$transcript_id[tt]))])
  }
  
  for (tt in c(1:length(finallist[,1])))
  {
    for (j in c(7:(dim(OC)[2]-1)))
    {
      finallist[tt,(j-3+length(c(7:(dim(OC)[2]-1))))]<-as.numeric(finallist[tt,j-4])*1000000000/(t[j-6]*as.numeric(finallist$length[tt]))
    }
    finallist$meanFPKM[tt]<-mean(as.numeric(finallist[tt,c((4+length(c(7:(dim(OC)[2]-1)))):(dim(OC)[2]-4+(dim(OC)[2]-5)))]))
  }
  
  
  
  length(finallist$gene_id[!duplicated(finallist$gene_id)]) #100
  finallist<-finallist[which(finallist$transcript_id %in% gtfcodingtrans$transcript_id),] 
  finallist<-finallist[which(finallist$transcript_id %in% gtfCDS$transcript_id),]
  length(finallist$gene_id[!duplicated(finallist$gene_id)]) #100
  finallist$ccdsid<-rep("0",dim(finallist)[1])
  for (j in c(1:dim(finallist)[1]))
  {
    tt=gtfCDS$ccdsid[which(gtfCDS$transcript_id==finallist$transcript_id[j])]
    finallist$ccdsid[j]<-tt[!duplicated(tt)]
  }
  finallist<-finallist[which(is.na(finallist$ccdsid)==FALSE),]
  length(finallist$gene_id[!duplicated(finallist$gene_id)]) #100
  max200gene<-finallist$gene_id[!duplicated(finallist$gene_id)]
  max200trans<-c()
  for (j in c(1:length(max200gene)))
  {
    ttt<-finallist[which(finallist$gene_id==max200gene[j]),]
    if (length(ttt$transcript_id[which(ttt$meanFPKM==max(ttt$meanFPKM))])==1)
    {
      max200trans[j]<-ttt$transcript_id[which(ttt$meanFPKM==max(ttt$meanFPKM))]
    }
    else
      if (length(ttt$transcript_id[which(ttt$meanFPKM==max(ttt$meanFPKM))])!=1)
      {
        max200trans[j]<-ttt$transcript_id[which(ttt$meanFPKM==max(ttt$meanFPKM))][1]
        print("warning")
        print(j)
      }
  }
  last<-data.frame(max200gene,max200trans)
  colnames(last)<-c("gene","transcript") #transcripts with highest FPKM
  last$protein<-rep("0",dim(last)[1])
  for (j in c(1:length(last[,1])))
  {
    t<-gtfCDS$protein_id[which(str_detect(gtfCDS$transcript_id,last$transcript[j])==TRUE)]
    last$protein[j]<-t[!duplicated(t)][1]
  }
  allgene<-last$gene
  alltranscript<-last$transcript
  allgeneproteinsequence<-list()
  for (tt in c(1:length(alltranscript)))
  {
    gtfCDSgene<-gtfCDS[which(gtfCDS$transcript_id==alltranscript[tt]),]
    proteinalternative<-gtfCDSgene$protein_id[!duplicated(gtfCDSgene$protein_id)]
    proteinsequence<-c()
    for (j in c(1:length(proteinalternative)))
    {
      gtfCDSgeneprotein<-gtfCDSgene[which(gtfCDSgene$protein_id==proteinalternative[j]),]
      CDS<-c()
      for (k in c(1:length(gtfCDSgeneprotein[,1])))
      {
        CDS[k]<-reverse(str_sub(fasta[chromosome[gtfCDSgeneprotein[k,1],1],2],gtfCDSgeneprotein[k,2],gtfCDSgeneprotein[k,3]),gtfCDSgeneprotein[k,5])
      }
      CDSseq<-pastevector(CDS,"") #CDS
      proteinsequence[j]<-CDSseq
    }
    names(proteinsequence)<-proteinalternative
    allgeneproteinsequence[[tt]]<-proteinsequence
    names(allgeneproteinsequence)[tt]<-alltranscript[tt]
  }
  identical(names(allgeneproteinsequence),last$transcript) #True
  last$CDSsequence<-unlist(allgeneproteinsequence)
  setwd("//home//haoyu//Desktop//PBMC//maxgene")
  write.csv(last,paste(celltype[i],"top200codinggenemaxtranscripts.csv",sep=""))
  dim(last)
  #finished
}