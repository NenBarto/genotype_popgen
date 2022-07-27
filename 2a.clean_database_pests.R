
library(taxonomizr)
library(Biostrings)
library(ggplot2)

projectname<-"18S"
infile=paste0("/share/ScratchGeneral/nenbar/projects/Anthony/scripts/",projectname,"/db/pests.11.fasta")

seqs<-readDNAStringSet(infile)
seqs<-seqs[order(names(seqs))]
#import data from table
annotation<-read.csv("db/MarinePestReferenceDatabase.csv")

df<-data.frame(ids=names(seqs))
df$ids<-gsub("~.*","",df$ids)
df$ids<-gsub("_.*","",df$ids)

merged<-merge(df,annotation,by.x="ids",by.y="Sample.ID",sort=F,all.x=T)


#plot widths
pdf("amplicon_width_pests.11.pdf",width=6,height=4)
qplot(width(seqs))+geom_histogram()+xlab("Amplicon width")
dev.off()

#convert all the ids to taxonomy
#seqsS<-seqs
#ids<-names(seqsS)
ids<-merged$Species.Name
idsClean<-gsub(" $","",ids)

#fix the names
#idsClean[1]<-"Dreissena rostriformis bugensis"
#idsClean[19:23]<-"Codium fragile"
#idsClean[117]<-"Pinctada albina"
#idsClean[124:125]<-"Caulerpa cylindracea"
#idsClean[126]<-"Pinctada albina"
#idsClean[c(139,141)]<-"Oestridae sp."

#fix 86


prepareDatabase('accessionTaxa.sql')


taxaId<-getId(idsClean,sqlFile="accessionTaxa.sql")

taxaId[1]=45949
taxaId[2]=427924
taxaId[29]=2589375
taxaId[35]=3133
taxaId[36]=3133
taxaId[37]=3133
taxaId[38]=3133
taxaId[39]=3133
taxaId[40]=3133
taxaId[42]=1289581
taxaId[149]=412131
taxaId[150]=2848398
taxaId[151]=2848398
taxaId[153]=1835367
taxaId[186]=315487
taxaId[191]=219692
taxaId[192]=219692
taxaId[193]=219692
taxaId[194]=2559386
taxaId[213]=1906867
taxaId[214]=123734
taxaId[216]=123734
taxaId[217]=123734
taxaId[218]=2790913
taxaId[219]=1295089



taxa<-getTaxonomy(taxaId,'accessionTaxa.sql',desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus",
         "species"))

#fix 153 by hand Charybdis yaldwyni
taxa[153,"species"]="Charybdis yaldwyni"
taxa[194,"species"]="Nectocarcinus antarcticus"
taxa[213,"species"]="Brachidontes sp."
taxa[c(214,216,217),"species"]="Oestridae sp."
taxa[218,"species"]="Xiphonectes rugosus"


tmp=taxa
row.names(tmp)<-1:dim(tmp)[1]


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


names(seqs)<-seqNames
writeXStringSet(seqs,"pests.11.taxed.fasta",format="fasta")








	