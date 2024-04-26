
### 24.04.26 compare binding locations of KDM5C ChIp-seq peaks and number of germline genes bound


# Visualize genomic locations using ChIPseeker
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
library(clusterProfiler)
library(ggplot2)
library(dplyr)



#get the parameter used for the promoter range
source("code/utilities/parameters.R")

#using WT no KO for both EpiLCs and PNCs
samples <- c("EpiLC", "PNC") #make sure order matches snakefile
#a list to save the plots in
annoplot <- list()
#a list of the peak gene names
peakENSEMBL <- list()


for (i in 1:length(samples)){
	#make a GRanges object from your .bed file of peak locations
	peak <- readPeakFile(snakemake@input[[i]])
	head(peak)

	#annotate where the peaks are located and plot
	peakAnno <- annotatePeak(peak, tssRegion=c(-promorange, promorange),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
	annoplot[[i]] <- plotAnnoBar(peakAnno, title = samples[i])


	#subset for just the promoter peaks
	annot_KDM5C <- as.data.frame(peakAnno@anno)
	annot_KDM5C.promo <- subset(annot_KDM5C, annotation == "Promoter")
	KDM5C.promo_genes <- annot_KDM5C.promo %>% distinct(ENSEMBL, .keep_all = TRUE) #keep only distinct annotations (in case 2 peaks are within 1000bp)
	#print(head(KDM5C.promo_genes))


	peakENSEMBL[[i]] <- unique(KDM5C.promo_genes$ENSEMBL)


}

library("ggplot2")
library("gridExtra")
library("ggpubr")
#ggsave(snakemake@output[[1]], grid.arrange(grobs = annoplot, nrow = 2), width = 6, height = 6)



#See which genes that kdm5c is bound to in epilcs and pncs overlap with DEGs
#put the name of the peak samples with their ensembls
names(peakENSEMBL) <- samples

#read in the DEGs

