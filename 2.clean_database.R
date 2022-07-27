
library(taxonomizr)
library(Biostrings)
library(ggplot2)
infile="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/split/trimmed_new/16S_trimmed.fasta"

seqs<-readDNAStringSet(infile)

#out of 167100 sequence ids there were 78698 sequences

#where same species and duplicated - collapse
#otherwise call ambiguous?

#plot widths
pdf("amplicon_width.pdf",width=6,height=4)
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
table(taxa[,2])
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

isInsect<-taxa[,3]=="Insecta"

names(seqsS)<-seqNames
writeXStringSet(seqsS,"insects.fasta",format="fasta")

writeXStringSet(seqsS[isInsect],"insectsOnly.fasta",format="fasta")


#








	
