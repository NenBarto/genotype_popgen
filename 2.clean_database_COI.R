
library(taxonomizr)
library(Biostrings)
library(ggplot2)

infile="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/COI/metaCOXI.fasta"

seqs<-readDNAStringSet(infile)

#compared to 16S:
#out of 92064 sequence ids there were 50564 unique sequences

#where same species and duplicated - collapse
#otherwise call ambiguous?

#plot widths
pdf("amplicon_width_metaCOXI_18S.pdf",width=6,height=4)
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
write.table(df,"db_tax_metaCOXI.11.txt",quote=F,sep="\t",row.names=F)

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

#isInsect<-taxa[,3]=="Insecta"

names(seqsS)<-seqNames
writeXStringSet(seqsS,"metaCOXI.11.taxed.fasta",format="fasta")









	