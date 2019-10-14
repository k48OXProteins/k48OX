setwd("~/Keira/")
edgeRout <- read.table("mek10_dn10.txt")
genes <- read.table("KH2013-UniqueNAME_def.txt",T,row.names = "GeneID")
gff <- read.table("6_10_18_annotated_bedtools_window_10000.txt")
peakIDs <- unique(gff[,4])
peakIDs <- head(unique(gff[,c(4,17)]))
geneIDs <- unique(gff[,17])
edgeRout$genes <- gff[,16][row.names(gff)]