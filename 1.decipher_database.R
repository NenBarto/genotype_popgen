
library(DECIPHER)
library(taxonomizr)
library(Biostrings)
library(ggplot2)


seqs_path<-"/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/out/inverts/trimmed_new/16S_trimmed.fasta"

seqs <- readDNAStringSet(seqs_path)
seqs <- RemoveGaps(seqs)
seqs <- OrientNucleotides(seqs)

#create tax ID file
seqnames<-gsub(";.*","",names(seqs))
idsClean<-gsub("\\..*","",seqnames)
prepareDatabase('accessionTaxa.sql')
taxaId<-accessionToTaxa(idsClean,sqlFile="accessionTaxa.sql",version='base')
taxonomy<-getTaxonomy(taxaId,'accessionTaxa.sql',desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus",
         "species"))

set.seed(123)
temp<-taxonomy

#clean from NAs
seqs<-seqs[!is.na(taxonomy[,"kingdom"])]

seqNames<-paste0(
	
	taxonomy[,"species"],
	" Root; ",taxonomy[,"kingdom"],
	"; ",taxonomy[,"phylum"],
	"; ",taxonomy[,"class"],
	"; ",taxonomy[,"order"],
	"; ",taxonomy[,"family"],
	"; ",taxonomy[,"genus"],
	"; ",taxonomy[,"species"])
seqNames<-seqNames[!is.na(taxonomy[,"kingdom"])]
names(seqs)<-seqNames

#taxa <- setNames(c("kingdom", "phylum", "class", "order", "family", "genus","species"),
#c("d__", "p__", "c__", "o__", "f__", "g__","s__"))

groups <- names(seqs)
groups <- gsub("(.*)(Root;)", "\\2", groups)
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups
length(u_groups) # number of groups

#53032 for 16S

#limit the classifier to 10 sequences
maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
	index <- which(groups==u_groups[i])
	keep <- sample(length(index),
	maxGroupSize)
	remove[index[-keep]] <- TRUE
}
sum(remove) # number of sequences eliminated

#removed 60466 sequences

maxIterations <- 3 # must be >= 1
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
	cat("Training iteration: ", i, "\n", sep="")
	# train the classifier
	trainingSet <- LearnTaxa(seqs[!remove],names(seqs)[!remove],taxid)
	# look for problem sequences
	probSeqs <- trainingSet$problemSequences$Index
	if (length(probSeqs)==0) {
	cat("No problem sequences remaining.\n")
	break
	} else if (length(probSeqs)==length(probSeqsPrev) &&
	all(probSeqsPrev==probSeqs)) {
	cat("Iterations converged.\n")
	break
	}

	if (i==maxIterations)
	break
	probSeqsPrev <- probSeqs
	
	# remove any problem sequences
	index <- which(!remove)[probSeqs]
	remove[index] <- TRUE # remove all problem sequences
	if (!allowGroupRemoval) {
		# replace any removed groups
		missing <- !(u_groups %in% groups[!remove])
		missing <- u_groups[missing]
		if (length(missing) > 0) {
			index <- index[groups[index] %in% missing]
			remove[index] <- FALSE # don't remove
		}
	}
}
sum(remove) # total number of sequences eliminated
length(probSeqs) # number of remaining problem sequences






ranks<-paste0(
	
	"d__",taxonomy[,"kingdom"],
	";p__",taxonomy[,"phylum"],
	";o__",taxonomy[,"order"],
	";f__",taxonomy[,"family"],
	";g__",taxonomy[,"genus"])

taxa <- setNames(c("domain", "phylum", "order", "family", "genus"),
c("d__", "p__", "o__", "f__", "g__"))



ranks <- strsplit(ranks, ";", fix=T)
count <- 1L
groups <- "Root"
index <- -1L
level <- 0L
rank <- "rootrank"
pBar <- txtProgressBar(style=3)

for (i in seq_along(ranks)) {
	for (j in seq_along(ranks[[i]])) {
		rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
		group <- substring(ranks[[i]][j], 4)
		w <- which(groups==group & rank==rank_level)
		if (length(w) > 0) {
			parent <- match(substring(ranks[[i]][j - 1], 4),
				groups)
			if (j==1 || any((parent - 1L)==index[w]))
				next # already included
		}
		count <- count + 1L
		groups <- c(groups, group)
		if (j==1) {
			index <- c(index, 0)
		} else {
			parent <- match(substring(ranks[[i]][j - 1], 4),groups)
		index <- c(index, parent - 1L)
		}
	level <- c(level, j)
	rank <- c(rank, taxa[j])
	}
	setTxtProgressBar(pBar, i/length(ranks))
}

groups <- gsub("^[ ]+", "", groups)
groups <- gsub("[ ]+$", "", groups)
taxid <- paste(0:(length(index) - 1L), groups, index, level, rank, sep="*")
head(taxid, n=10)
taxidFile="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/16S_trimmednew_taxid.txt"
writeLines(taxid,
con=taxidFile)

taxid <- read.table(taxidFile,
	header=FALSE,
	col.names=c('Index', 'Name', 'Parent', 'Level', 'Rank'),
	sep="*", # asterisks delimited
	quote="", # preserve quotes
	stringsAsFactors=FALSE)









































maxIterations <- 3 
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() 

for (i in seq_len(maxIterations)) {
	cat("Training iteration: ", i, "\n", sep="")
	# train the classifier
	trainingSet <- LearnTaxa(seqs,
	names(seqs),
	taxid)
	# look for problem sequences
	probSeqs <- trainingSet$problemSequences$Index
	if (length(probSeqs)==0) {
	cat("No problem sequences remaining.\n")
	break
	} else if (length(probSeqs)==length(probSeqsPrev) &&
	all(probSeqsPrev==probSeqs)) {
	cat("Iterations converged.\n")
	break
}



indir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/size10_0.8_80leftover/"
infasta=paste0(indir,"all.otus.fasta")

seqs<-readDNAStringSet(infasta)
seqnames<-gsub(";.*","",names(seqs))
names(seqs)<-seqnames
#there are 24678 OTUs
#now get those that are unassigned

unassigned<-read.table(paste0(indir,"counts_with_tax.txt"),sep="\t",header=T)

unassignedOTUs<-unassigned$OTU_ID[unassigned$TAX=="Unassigned"]
#21173 out of 24502 are unassigned previously
#14955 out of 24678 are unassigned now

seqsShort<-seqs[names(seqs) %in% unassignedOTUs]

#run blast
system("conda activate enviroDNA")
system(paste0("blastn -query ",indir,"/all.otus.fasta -db ../../annotation/blast/16S_metazoa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qcovs evalue bitscore\" -out all.otus-16S_metazoa.tab"))
system(paste0("mv all.otus-16S_metazoa.tab ",indir))

#now load in the results of the blast hits and filter them
blastDF<-read.table(paste0(indir,"all.otus-16S_metazoa.tab"))
blastDF$id<-gsub(";.*","",blastDF$V1)
blastDF$gb<-gsub("\\|$","",blastDF$V2)
blastDF$gb<-gsub(".*\\|","",blastDF$gb)

blastDFS<-blastDF[blastDF$V3<=80,]
#take the best hit per id, and only the first one 
#this could be fixed by taking multiple and then taking unique taxa
#fix later

blastDFSL<-split(blastDFS,blastDFS$id)
blastGB<-lapply(blastDFSL,function(x){x[,c("id","gb")][x$V3==max(x$V3),][1,]})
blastGBdf<-do.call("rbind",blastGB)

#now convert the ids
ids<-names(seqsShort)
df<-data.frame(otusOld=ids)
merged<-merge(df,blastGBdf,by.x="otusOld",by.y="id",sort = F)
#convert all the ids to taxonomy
seqsS<-seqsShort[names(seqsShort) %in% merged$otusOld]
names(seqsS)=merged$gb

