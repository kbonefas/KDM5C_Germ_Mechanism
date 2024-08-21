# 24.08.19 - WGBS volcano plots

#read in the KDM5C bound genes
KDM5C_binding <- read.csv(snakemake@input[[1]], sep = ",")

print(head(KDM5C_binding))

#get the gene coordinate information
source("code/utilities/GeneTSSandTES.R")
germTSS <- geneTSS_df(KDM5C_binding$ENSEMBL, 500)
print(head(germTSS))


samples <- c("all germ promoters", "germ CGI promoters")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(ggplot2)
library(dplyr)

##### #make the dataframes for the volcano plot
## TSS
TSSdf <- read.csv(snakemake@input[[2]], sep = ",")

#rename chr to seqnames to match granges for merging later
names(TSSdf)[names(TSSdf) == "chr"]<- "seqnames"
print(paste("rows in all germ -", nrow(TSSdf)))
print(head(TSSdf))
	
#get the germline genes for each coordinate
plotTSS <- merge(germTSS, TSSdf)
plotTSS <- merge(KDM5C_binding, plotTSS, by = "ENSEMBL")
#remove any duplicated loci
plotTSS <- plotTSS[!duplicated(plotTSS[6:7]),]
print(paste("rows in merged all germ -", nrow(plotTSS)))
print(head(plotTSS))

names(plotTSS)[names(plotTSS) == "gene_name"] <- "SYMBOL"



write.table(plotTSS, snakemake@output[[1]], sep = ",", row.names = FALSE)






#####
## CGIs
CGIdf <- read.csv(snakemake@input[[3]], sep = ",")

#rename chr to seqnames to match granges for merging later
names(CGIdf)[names(CGIdf) == "chr"]<- "seqnames"
print(paste("rows in CGI -", nrow(CGIdf)))
print(head(CGIdf))
	
#get the germline genes for each coordinate
	#make granges object with bed coordinates
peak <- makeGRangesFromDataFrame(CGIdf)
print(head(peak))

#annotate locations
peakAnno <- annotatePeak(peak, tssRegion=c(-500, 500),
                    TxDb=txdb, annoDb="org.Mm.eg.db")
	
annot <- as.data.frame(peakAnno@anno)
annot_promo <- subset(annot, annotation == "Promoter")
#annot_promo_genes <- annot_promo %>% distinct(ENSEMBL, .keep_all = TRUE)

print(paste("rows in annotated CGIs -", nrow(annot_promo)))
print(head(annot_promo))

plotCGI <- merge(annot_promo, CGIdf)
plotCGI <- subset(plotCGI, select = c(seqnames, start, end, ENSEMBL, SYMBOL))
plotCGI <- merge(KDM5C_binding, plotCGI)
print(paste("rows in CGI merged -", nrow(plotCGI)))
print(head(plotCGI))

write.table(plotCGI, snakemake@output[[2]], sep = ",",  row.names = FALSE)


