library(tidyverse)
library(gtools)

args=commandArgs(trailingOnly=TRUE)

if (length(args)!=2)
{
  stop("2 args [input-file] [output-file] must be supplied", call.=FALSE)
}

d1<-read.delim(args[1], header=TRUE)
d1<-d1[,mixedsort(colnames(d1))] #from gtools

#d1<-read.delim("/Users/svs/Downloads/ReanalyseFish/vsearch/id100-minfreq8-sintax100-withcarp/counts_with_tax.txt", header=TRUE)

#select syntax comes from tidyverse package, - means exclude, %>% means pipe
sum_tab = 
  d1 %>% 
  select(-OTU_ID) %>% 
  group_by(TAX) %>% 
  #summarise_all(funs(sum))
  summarise_all(list(~sum(.)))
#View(sum_tab)

write.table(sum_tab, file=args[2], sep=";", row.names=FALSE)