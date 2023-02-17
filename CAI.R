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

Generatecodon<-function(group=TRUE) #subgroup Ser,Leu,Arg or not
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
} #generate codon table

RSCU<-function(highexpression,group=TRUE) #Input example: highexpression=c(pastevector(rownames(codonaas),""),"AGGTTTTATTAT"), from ATG to TAG/TAA/TGA
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


CAI<-function(sequence,codontable,adj=TRUE) #Input example: sequence="AGGTTTATG"
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



#RSCU


setwd("~//CAImaxgeneCDS")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImaxgeneCDS",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  setwd("~\\CAI")
  write.csv(codontable,paste(paste(celltypelist[i],"top200codinggeneRSCU",sep=""),"csv",sep="."))
}

setwd("~//CAImaxgeneCDS//PBMC//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//PBMC//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  setwd("~\\CAI")
  write.csv(codontable,paste(paste(paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),"top200codinggeneRSCU",sep=""),"csv",sep="."))
}  

setwd("~//CAImaxgeneCDS//lung//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//lung//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  setwd("~\\CAI")
  write.csv(codontable,paste(paste(paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),"top200codinggeneRSCU",sep=""),"csv",sep="."))
}


setwd("~//CAImxgeneCDS(afterviralinfection)")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  setwd("~\\CAI")
  write.csv(codontable,paste(paste(celltypelist[i],"top200codinggeneRSCU",sep=""),"csv",sep="."))
}



setwd("~//CAImxgeneCDS(afterviralinfection)new")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)new",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  setwd("~\\CAI")
  write.csv(codontable,paste(paste(celltypelist[i],"top200codinggeneRSCU",sep=""),"csv",sep="."))
}




###HIV

setwd("~//CAI//HIV")
ALLHIV<-list.files()
ALLHIV<-ALLHIV[str_detect(ALLHIV,".csv")]
ALLHIV<-ALLHIV[str_detect(ALLHIV,"HIVvariantsequence")==FALSE]
setwd("~//CAImaxgeneCDS")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
print(celltypelist[i])
setwd(paste("~//CAImaxgeneCDS",celltypelist[i],sep="//"))
TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
codontable<-RSCU(TOP200,group=TRUE)
for (j in c(1:length(ALLHIV)))
{
setwd("~//CAI//HIVsequence")
HIVCDS<-read.csv(ALLHIV[j])
CAIlist<-c()
CAIlistadj<-c()
for (t in c(1:dim(HIVCDS)[1]))
{
if (HIVCDS$normallist[t]!="function")
{
  CAIlist[t]<-"inapplicability"
  CAIlistadj[t]<-"inapplicability"
}
if (HIVCDS$normallist[t]=="function")
{
  CAIlist[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=FALSE)
  CAIlistadj[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=TRUE)
}
}
HIVCDS$CAI<-CAIlist
HIVCDS$CAIadj<-CAIlistadj
setwd("~//CAI//HIVsequenceCDS")
dir.create(celltypelist[i])
setwd(paste("~//CAI//HIVsequenceCDS",celltypelist[i],sep="//"))
write.csv(HIVCDS[,-1],paste(paste(ALLHIV[j],celltypelist[i],sep=" in "),"csv",sep="."))
setwd("~//CAI//HIVsequence")
}
}


setwd("~//CAImaxgeneCDS//PBMC//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//PBMC//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLHIV)))
  {
    setwd("~//CAI//HIVsequence")
    HIVCDS<-read.csv(ALLHIV[j])
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(HIVCDS)[1]))
    {
      if (HIVCDS$normallist[t]!="function")
      {
        CAIlist[t]<-"inapplicability"
        CAIlistadj[t]<-"inapplicability"
      }
      if (HIVCDS$normallist[t]=="function")
      {
        CAIlist[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=TRUE)
      }
    }
    HIVCDS$CAI<-CAIlist
    HIVCDS$CAIadj<-CAIlistadj
    setwd("~//CAI//HIVsequenceCDS")
    dir.create(paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//HIVsequenceCDS",paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    write.csv(HIVCDS[,-1],paste(paste(ALLHIV[j],paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//HIVsequence")
  }
}

setwd("~//CAImaxgeneCDS//lung//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//lung//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLHIV)))
  {
    setwd("~//CAI//HIVsequence")
    HIVCDS<-read.csv(ALLHIV[j])
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(HIVCDS)[1]))
    {
      if (HIVCDS$normallist[t]!="function")
      {
        CAIlist[t]<-"inapplicability"
        CAIlistadj[t]<-"inapplicability"
      }
      if (HIVCDS$normallist[t]=="function")
      {
        CAIlist[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=TRUE)
      }
    }
    HIVCDS$CAI<-CAIlist
    HIVCDS$CAIadj<-CAIlistadj
    setwd("~//CAI//HIVsequenceCDS")
    dir.create(paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//HIVsequenceCDS",paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    write.csv(HIVCDS[,-1],paste(paste(ALLHIV[j],paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//HIVsequence")
  }
}

setwd("~//CAI//HIVsequence")
ALLHIV<-list.files()
ALLHIV<-ALLHIV[str_detect(ALLHIV,".csv")]
ALLHIV<-ALLHIV[str_detect(ALLHIV,"HIVvariantsequence")==FALSE]
setwd("~//CAImxgeneCDS(afterviralinfection)")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLHIV)))
  {
    setwd("~//CAI//HIVsequence")
    HIVCDS<-read.csv(ALLHIV[j])
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(HIVCDS)[1]))
    {
      if (HIVCDS$normallist[t]!="function")
      {
        CAIlist[t]<-"inapplicability"
        CAIlistadj[t]<-"inapplicability"
      }
      if (HIVCDS$normallist[t]=="function")
      {
        CAIlist[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=TRUE)
      }
    }
    HIVCDS$CAI<-CAIlist
    HIVCDS$CAIadj<-CAIlistadj
    setwd("~//CAI//HIVsequenceCDS(afterviralinfection)")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//HIVsequenceCDS(afterviralinfection)",celltypelist[i],sep="//"))
    write.csv(HIVCDS[,-1],paste(paste(ALLHIV[j],celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//HIVsequence")
  }
}

setwd("~//CAI//HIVsequence")
ALLHIV<-list.files()
ALLHIV<-ALLHIV[str_detect(ALLHIV,".csv")]
ALLHIV<-ALLHIV[str_detect(ALLHIV,"HIVvariantsequence")==FALSE]
setwd("~//CAImxgeneCDS(afterviralinfection)new")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)new",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLHIV)))
  {
    setwd("~//CAI//HIVsequence")
    HIVCDS<-read.csv(ALLHIV[j])
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(HIVCDS)[1]))
    {
      if (HIVCDS$normallist[t]!="function")
      {
        CAIlist[t]<-"inapplicability"
        CAIlistadj[t]<-"inapplicability"
      }
      if (HIVCDS$normallist[t]=="function")
      {
        CAIlist[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(HIVCDS$sequencelist[t],codontable,adj=TRUE)
      }
    }
    HIVCDS$CAI<-CAIlist
    HIVCDS$CAIadj<-CAIlistadj
    setwd("~//CAI//HIVsequenceCDS(afterviralinfection)Ribo")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//HIVsequenceCDS(afterviralinfection)Ribo",celltypelist[i],sep="//"))
    write.csv(HIVCDS[,-1],paste(paste(ALLHIV[j],celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//HIVsequence")
  }
}














###SARSCOV2
library(Biostrings)
library(ggplot2)
library(stringr)
coronagenes<-c("ORF1ab","ORF1a","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
setwd("~//CAI//SARSCOV2sequence")
ALLCORONA<-list.files()
ALLCORONA<-ALLCORONA[str_detect(ALLCORONA,".fasta")]
setwd("~//CAImaxgeneCDS")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImaxgeneCDS",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//SARSCOV2sequence")
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//SARSCOV2sequenceCDS")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//SARSCOV2sequenceCDS",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    #write.csv(CORONACDS[,-1],paste(paste(paste(str_split(str_split(ALLCORONA[j],"\\(")[[1]][2],"\\)")[[1]][1],str_split(names(fasta_input)[1]," ")[[1]][1],sep="-"),celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//SARSCOV2sequence")
  }
}


setwd("~//CAImaxgeneCDS//PBMC//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//PBMC//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//SARSCOV2sequence")
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
        CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//SARSCOV2sequenceCDS")
    dir.create(paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//SARSCOV2sequenceCDS",paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    #write.csv(CORONACDS[,-1],paste(paste(paste(str_split(str_split(ALLCORONA[j],"\\(")[[1]][2],"\\)")[[1]][1],str_split(names(fasta_input)[1]," ")[[1]][1],sep="-"),paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//SARSCOV2sequence")
  }
}

setwd("~//CAImaxgeneCDS//lung//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//lung//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//SARSCOV2sequence")
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
        CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
        CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//SARSCOV2sequenceCDS")
    dir.create(paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//SARSCOV2sequenceCDS",paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    #write.csv(CORONACDS[,-1],paste(paste(paste(str_split(str_split(ALLCORONA[j],"\\(")[[1]][2],"\\)")[[1]][1],str_split(names(fasta_input)[1]," ")[[1]][1],sep="-"),paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//SARSCOV2sequence")
  }
}

setwd("~//CAI//SARSCOV2sequence")
ALLCORONA<-list.files()
ALLCORONA<-ALLCORONA[str_detect(ALLCORONA,".fasta")]
setwd("~//CAImxgeneCDS(afterviralinfection)")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//SARSCOV2sequence")
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//SARSCOV2sequenceCDS(afterviralinfection)")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//SARSCOV2sequenceCDS(afterviralinfection)",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    #write.csv(CORONACDS[,-1],paste(paste(paste(str_split(str_split(ALLCORONA[j],"\\(")[[1]][2],"\\)")[[1]][1],str_split(names(fasta_input)[1]," ")[[1]][1],sep="-"),celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//SARSCOV2sequence")
  }
}

setwd("~//CAI//SARSCOV2sequence")
ALLCORONA<-list.files()
ALLCORONA<-ALLCORONA[str_detect(ALLCORONA,".fasta")]
setwd("~//CAImxgeneCDS(afterviralinfection)new")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)new",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//SARSCOV2sequence")
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//SARSCOV2sequenceCDS(afterviralinfection)Ribo")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//SARSCOV2sequenceCDS(afterviralinfection)Ribo",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    #write.csv(CORONACDS[,-1],paste(paste(paste(str_split(str_split(ALLCORONA[j],"\\(")[[1]][2],"\\)")[[1]][1],str_split(names(fasta_input)[1]," ")[[1]][1],sep="-"),celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//SARSCOV2sequence")
  }
}














###othercoronavirus
library(Biostrings)
library(ggplot2)
library(stringr)
setwd("~//CAI//othercoronavirussequence")
coronageneslist<-list()
coronageneslist[[1]]<-c("ORF1ab","ORF1a","S","ORF4a","ORF4b","E","M","N")
coronageneslist[[2]]<-c("ORF1ab","ORF1a","HE","S","ORF4","E","M","N","N2")
coronageneslist[[3]]<-c("ORF1ab","ORF1a","S","ORF3","E","M","N")
coronageneslist[[4]]<-c("ORF1ab","ORF1a","NS2","HE","S","NS12.9","E","M","N","N2")
coronageneslist[[5]]<-c("ORF1ab","ORF1a","S","ORF3","ORF4a","ORF4b","ORF5","E","M","N","ORF8b")
coronageneslist[[6]]<-c("ORF1ab","ORF1a","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
coronageneslist[[7]]<-c("ORF1ab","ORF1a","S","ORF3a","ORF3b","E","M","ORF6","ORF7a","ORF7b","ORF8a","ORF8b","N","ORF9b","ORF9c")
ALLCORONA<-list.files()
ALLCORONA<-ALLCORONA[str_detect(ALLCORONA,".fasta")]
setwd("~//CAImaxgeneCDS")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImaxgeneCDS",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//othercoronavirussequence")
    coronagenes<-coronageneslist[[j]]
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//othercoronavirussequenceCDS")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//othercoronavirussequenceCDS",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//othercoronavirussequence")
  }
}


setwd("~//CAImaxgeneCDS//PBMC//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//PBMC//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//othercoronavirussequence")
    coronagenes<-coronageneslist[[j]]
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//othercoronavirussequenceCDS")
    dir.create(paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//othercoronavirussequenceCDS",paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],paste("PBMC",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//othercoronavirussequence")
  }
}

setwd("~//CAImaxgeneCDS//lung//maxgene")
celltypelist<-list.files()
for(i in c(1:length(celltypelist)))
{ 
  setwd("~//CAImaxgeneCDS//lung//maxgene")
  print(celltypelist[i])
  TOP200<-read.csv(celltypelist[i])[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//othercoronavirussequence")
    coronagenes<-coronageneslist[[j]]
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//othercoronavirussequenceCDS")
    dir.create(paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""))
    setwd(paste("~//CAI//othercoronavirussequenceCDS",paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],paste("lung",str_split(celltypelist[i],"top200")[[1]][1],sep=""),sep=" in "),"csv",sep="."))
    setwd("~//CAI//othercoronavirussequence")
  }
}

setwd("~//CAImxgeneCDS(afterviralinfection)")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//othercoronavirussequence")
    coronagenes<-coronageneslist[[j]]
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//othercoronavirussequenceCDS(afterviralinfection)")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//othercoronavirussequenceCDS(afterviralinfection)",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//othercoronavirussequence")
  }
}

setwd("~//CAImxgeneCDS(afterviralinfection)new")
celltypelist<-list.files()
celltypelist<-celltypelist[which(celltypelist %in% c("PBMC","lung")==FALSE)]
for(i in c(1:length(celltypelist)))
{ 
  print(celltypelist[i])
  setwd(paste("~//CAImxgeneCDS(afterviralinfection)new",celltypelist[i],sep="//"))
  TOP200<-read.csv("top200codinggenemaxtranscripts.csv")[,5]
  codontable<-RSCU(TOP200,group=TRUE)
  for (j in c(1:length(ALLCORONA)))
  {
    setwd("~//CAI//othercoronavirussequence")
    coronagenes<-coronageneslist[[j]]
    fasta_input<-readBStringSet(ALLCORONA[j], format='fasta')
    CORONACDS<-data.frame(coronagenes,width(fasta_input),as.data.frame(fasta_input)[,1])
    colnames(CORONACDS)<-c("genelist","CDSlength","sequencelist")
    CAIlist<-c()
    CAIlistadj<-c()
    for (t in c(1:dim(CORONACDS)[1]))
    {
      CAIlist[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=FALSE)
      CAIlistadj[t]<-CAI(CORONACDS$sequencelist[t],codontable,adj=TRUE)
    }
    CORONACDS$CAI<-CAIlist
    CORONACDS$CAIadj<-CAIlistadj
    setwd("~//CAI//othercoronavirussequenceCDS(afterviralinfection)Ribo")
    dir.create(celltypelist[i])
    setwd(paste("~//CAI//othercoronavirussequenceCDS(afterviralinfection)Ribo",celltypelist[i],sep="//"))
    write.csv(CORONACDS[,-1],paste(paste(str_split(ALLCORONA[j],"\\.")[[1]][1],celltypelist[i],sep=" in "),"csv",sep="."))
    setwd("~//CAI//othercoronavirussequence")
  }
}

















