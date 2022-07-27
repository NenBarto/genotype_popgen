
#this script tests bam files
library(GenomicAlignments)

filename="../../results/210813_miseq_AR.bwa/7111_S29.sorted.bam"
gr<-readGAlignmentPairs(filename)

first<-first(gr)
last<-last(gr)

df<-data.frame(fs=start(first),fe=end(first),ls=start(last),le=end(last))
all=apply(df,1,function(x){paste(x,collapse="_")})
length(unique(all))