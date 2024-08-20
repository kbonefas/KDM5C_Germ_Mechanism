# 24.08.19 - WGBS volcano plots

#read in the KDM5C bound genes
KDM5C_binding <- read.csv(snakemake@input[[1]], sep = ",")

print(head(KDM5C_binding))


samples <- c("all germ promoters", "germ CGI promoters")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(ggplot2)
library(dplyr)

volcano_WGBS <- function(n){
	df <- read.csv(snakemake@input[[n+1]], sep = ",")

	#rename chr to seqnames to match granges for merging later
	names(df)[names(df) == "chr"]<- "seqnames"
	print(paste("rows in", samples[n], "-", nrow(df)))
	print(head(df))
	

	#get the germline genes for each coordinate
		#make granges object with bed coordinates
	peak <- makeGRangesFromDataFrame(df)
	print(head(peak))

		#annotate locations
	peakAnno <- annotatePeak(peak, tssRegion=c(-500, 500),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
	
	annot <- as.data.frame(peakAnno@anno)
	annot_promo <- subset(annot, annotation == "Promoter")
	#annot_promo_genes <- annot_promo %>% distinct(ENSEMBL, .keep_all = TRUE)

	print(paste("rows in annotated", samples[n], "-", nrow(annot_promo)))
	print(head(annot_promo))

	merge <- merge(annot_promo, df)
	print(paste("rows in merged", samples[n], "-", nrow(merge)))
	print(head(merge))

	#its working for the CGIs but all of the promoters there are extras. maybe its because multiple genes have TSSs in that location
		#Do I necessarily need the gene location for the volcano plot? no. Just need a unique name for each region

}


for (i in 1:length(samples)){
	volcano_WGBS(i)
}


