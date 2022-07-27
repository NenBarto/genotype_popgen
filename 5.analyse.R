library(vcfR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(ggplot2)
library(pheatmap)
library(ggrepel)

projectdir="/share/ScratchGeneral/nenbar/projects/Anthony"
scriptsdir=paste0(projectdir,"/scripts")
tabledir=paste0(projectdir,"/tables")
system(paste0("mkdir -p ",tabledir))
figdir=paste0(projectdir,"/figures")
system(paste0("mkdir -p ",figdir))

#import vcf file
vcf.fn<-paste0(projectdir,"/results/EBB.vcf/mitochondrial.vcf.gz")
snpgdsVCF2GDS(vcf.fn, "EBB.gds")
genofile <- snpgdsOpen("EBB.gds")

#get genotypes and sample IDs
g <- read.gdsn(index.gdsn(genofile, "genotype"))
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))


snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only=F)
snpset.id <- unlist(unname(snpset))
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2,autosome.only=F)

#clean genotypes
gS<-as.data.frame(g)
row.names(gS)<-sample.id



pop_code <- gsub(".*\\.","",sample.id)
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)

pdf(paste0(figdir,"/PCA_EBB_mitochondrial.pdf"),height=10,width=10)
p<-ggplot(tab,aes(EV1,EV2,color=pop))
#p<-p+geom_jitter(width = 0.05, height = 0.05)
p<-p+geom_point()
p<-p+geom_text_repel(label=sample.id,max.overlaps=20)
p
dev.off()

pdf(paste0(figdir,"/heatmap_EBB_mitochondrial.pdf"),height=8,width=12)
pheatmap(gS)
dev.off()




