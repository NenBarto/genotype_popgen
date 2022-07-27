#this script takes 
library(rtracklayer)
data=read.table("../../annotation/EBB/EBB.scaffolds.fa.fai")
splitList<-cut(1:dim(data)[1],100)
gr<-GRanges(seqnames=data$V1,IRanges(1,data$V2))
#in the future do this to allow for equal size of vcf files and fewer failed processes
gr<-gr[sample(1:length(gr))]
grL<-split(gr,splitList)

grLE<-endoapply(grL,sort)
for(i in 1:length(grLE)){
  export.bed(grLE[[i]],paste0("EBB_chrs_",i,".bed"))
}