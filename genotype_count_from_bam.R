
library(Rsamtools)
library(vcfR)
inDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/210813_miseq_AR.bwa"
vcf_file<-"/share/ScratchGeneral/nenbar/projects/Anthony/results/finalVCF/raw.g5mac3dp10.recode.vcf"
vcf_file<-"/share/ScratchGeneral/nenbar/projects/Anthony/results/finalVCF/DP10g95maf05.recode.vcf"
vcf <- read.vcfR( vcf_file, verbose = FALSE )

positions<-as.numeric(vcf@fix[,2])

results<-list()
files<-list.files(inDir,pattern=".bam",full.names=T)
files<-files[!grepl("bai",files)]
for(file in files){
	sampleName<-gsub("_.*","",basename(file))
	cat(".")
	aln <- scanBam(file)
	seqs<-aln[[1]]$seq
	seqs<-seqs[width(seqs)==282]
	short<-sapply(seqs,function(x){as.character(x[positions])})
	results[[sampleName]]<-short
}

lengths<-sapply(results,length)
results=results[lengths>0]
uniqueList<-sapply(results,function(x){length(unique(x))})
sort(uniqueList,decreasing=T)
tableList<-lapply(results,table)

tableMean<-lapply(tableList,function(x){x/sum(x)})
tableMeanCounts<-lapply(tableMean,function(x){length(x[x>0.01])})
table(unlist(tableMeanCounts))


tableMeanCounts<-lapply(tableMean,function(x){length(x[x>0.03])})
table(unlist(tableMeanCounts))

tableMeanCounts<-lapply(tableMean,function(x){length(x[x>0.05])})
table(unlist(tableMeanCounts))

#Questions

#Population wise what are the most common genotypes and how are they spread across

#1. create the table with OTUs
#-take sequences 
#-sort by most common
#-convert to ID 
#2. measure how many time each OTU per population
#
#3. vcfR
results<-results[!grepl("negative",names(results))]
OTUs<-unlist(results)
df<-data.frame(OTU=OTUs,sampleName=rep(names(results),times=sapply(results,length)))
df$count<-1
dfAgg<-aggregate(df$count,list(OTU=df$OTU,sampleName=df$sampleName),sum)
colnames(dfAgg)<-c("OTU","sampleName","value")




aggregate(x$val1, list(id11 = x$id1, id2= x$id2), mean)
dfL<-split(df$sampleName,df$OTU)

dfTable<-lapply(dfL,table)

do.call("rbind",dfTable)

for(sampleName in unique(df$sampleName)){



}















