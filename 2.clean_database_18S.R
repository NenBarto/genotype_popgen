
library(taxonomizr)
library(Biostrings)
library(ggplot2)

#input from insilico PCR, 11 mismatches allowed
infile="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/18S/zhan/zhan.11.fasta"
seqs<-readDNAStringSet(infile)

#plot widths
pdf("amplicon_width_zhan11_18S.pdf",width=6,height=4)
qplot(width(seqs))+geom_histogram()+xlab("Amplicon width")
dev.off()

#convert all the ids to taxonomy
seqsS<-seqs
ids<-names(seqsS)
idsClean<-gsub("\\..*","",ids)
prepareDatabase('accessionTaxa.sql')


taxaId<-accessionToTaxa(idsClean,sqlFile="accessionTaxa.sql",version='base')
taxa<-getTaxonomy(taxaId,'accessionTaxa.sql',desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus",
         "species"))
df<-as.data.frame(table(taxa[,2]))
colnames(df)<-c("phylum","count")
write.table(df,"db_tax_zhan.11.txt",quote=F,sep="\t",row.names=F)

seqNames<-paste0(
	taxa[,"species"],"-",idsClean,
	";tax=k:",taxa[,"kingdom"],
	",p:",taxa[,"phylum"],
	",c:",taxa[,"class"],
	",o:",taxa[,"order"],
	",f:",taxa[,"family"],
	",g:",taxa[,"genus"],
	",s:",taxa[,"species"],
	"-",idsClean)
seqNames<-gsub(" ","-",seqNames)

names(seqsS)<-seqNames
writeXStringSet(seqsS,"zhan.11.taxed.fasta",format="fasta")









	
