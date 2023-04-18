# CAIvirus

upstream.sh: download datasets and get expression matrix at gene and transcripts level

gettop200maxgenelist.R:  get top200 highly-expression genes in each dataset

gettop200CDSlist.R:  get CDS sequence of transcripts with highest FPKM of top200 highly-expression genes from gettop200maxgenelist.R

scRNAseq-lung/PBMC.R: analyze two scRNA-seq datasets, cluster and annotate cells

CAI.R: calculate CAI of viral ORFs in each dataset

CAIandTEregression.R: calculate the relationship between CAI and translational efficiency

HIVanalysis: downstream analysis in HIV-1 part
