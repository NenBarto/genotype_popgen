#I first downloaded the sequences from NCBI with R package rentrez
#The sequence file, each of 100.000 fasta sequences were then split
#into chunks of 5000 sequences. 
#Each chunk was then processed with Dave's script


#download from
library(R.utils)
library(rentrez)
args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["param"]])){param = as.numeric(args$param)}
if (!is.null(args[["startnew"]])){startnew = as.numeric(args$startnew)}

cat(param)

metazoa_COI <- entrez_search(db="nuccore", term="(COI) AND mitochondrial) AND 1 : 10000[Sequence Length] AND animals[filter] ", use_history=TRUE,api_key ="60b6e6ed96096ce6ea0f8c049b625ffb6a09")

outdir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/COI/ncbi"

i=param
for( seq_start in seq(i*100000+startnew+1,c(i+1)*100000,250)){
	cat(seq_start)
	cat("\n")
    recs <- entrez_fetch(db="nuccore", web_history=metazoa_COI$web_history,
                         rettype="fasta", retmax=250, retstart=seq_start,api_key ="60b6e6ed96096ce6ea0f8c049b625ffb6a09")
    cat(recs, file=paste0(outdir,"/COI_",i,".fasta"), append=TRUE)
    cat(seq_start+249, "sequences downloaded\r")
}


