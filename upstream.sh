#bulk RNA-seq as example

#!/bin/bash

cd ~/target
R
a<-readLines("dataset")[1:6]
a1<-strsplit(a[1],",|;")[[1]]
a1[1]<-strsplit(a1[1]," ")[[1]][length(strsplit(a1[1]," ")[[1]])]
write.csv(data.frame(c(a1)),'pbmclist1.txt',row.names=FALSE)
q()
no
cat pbmclist1.txt | tail -n +2 > pbmclist.txt
rm pbmclist1.txt
sed -i 's/"/ /g' pbmclist.txt 

CC=`cat pbmclist.txt`

echo $CC
for i in $CC
do
prefetch $i
fastq-dump --split-3 --gzip ~/target/$i/$i.sra -O ~/target/$i/
cd ~/target/$i/
java -jar ~/Desktop/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 ~/target/$i/${i}_1.fastq.gz ~/target/$i/${i}_2.fastq.gz ~/target/$i/${i}_1.paired.fastq.gz ~/target/$i/${i}_2.paired.fastq.gz ~/target/$i/${i}_1.unpaired.fastq.gz ~/target/$i/${i}_2.unpaired.fastq.gz ILLUMINACLIP:/home/haoyu/Desktop/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
fastqc -t 6 -o ~/target/$i ~/target/$i/${i}_1.paired.fastq.gz
fastqc -t 6 -o ~/target/$i ~/target/$i/${i}_2.paired.fastq.gz
hisat2 -p 16 -x ~/Desktop/grch38/grch38/genome -1 ~/target/$i/${i}_1.paired.fastq.gz -2 ~/target/$i/${i}_1.unpaired.fastq.gz -S ~/target/$i/$i.sam --summary-file ~/target/$i/$i.hisat2.summary
samtools sort -@ 6 -o ~/target/$i/$i.bam ~/target/$i/$i.sam
rm ~/target/$i/$i.sam
featureCounts -p -T 8 -t gene -g gene_name -a ~/Desktop/gencode.v39.annotation.gtf  -o ~/target/$i.out ~/target/$i/$i.bam
featureCounts -p -T 8 -t exon -g gene_name -a ~/Desktop/gencode.v39.annotation.gtf  -o ~/target/$i.exon.out ~/target/$i/$i.bam
featureCounts -f -O -p -T 8 -F GTF -a  ~/Desktop/gencode.v39.annotation_dexseq.gtf  -o  ~/target/AS/$i.AS.out ~/target/$i/$i.bam
done

exit 0


