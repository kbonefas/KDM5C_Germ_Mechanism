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
# print(paste("rows in all germ -", nrow(TSSdf)))
# print(head(TSSdf))
	
#get the germline genes for each coordinate
plotTSS <- merge(germTSS, TSSdf)
plotTSS <- merge(KDM5C_binding, plotTSS, by = "ENSEMBL")
#remove any duplicated loci
plotTSS <- plotTSS[!duplicated(plotTSS[6:7]),]
#if there are duplicate gene names, 
print(paste("rows in merged all germ -", nrow(plotTSS)))
print(head(plotTSS))

#names(plotTSS)[names(plotTSS) == "gene_name"] <- "SYMBOL"


write.table(plotTSS, snakemake@output[[1]], sep = ",", row.names = FALSE)






#####
## CGIs
CGIdf <- read.csv(snakemake@input[[3]], sep = ",")

#rename chr to seqnames to match granges for merging later
names(CGIdf)[names(CGIdf) == "chr"]<- "seqnames"
# print(paste("rows in CGI -", nrow(CGIdf)))
# print(head(CGIdf))
	
#get the germline genes for each coordinate
	#make granges object with bed coordinates
peak <- makeGRangesFromDataFrame(CGIdf)
# print(head(peak))

#annotate locations
peakAnno <- annotatePeak(peak, tssRegion=c(-500, 500),
                    TxDb=txdb, annoDb="org.Mm.eg.db")
	
annot <- as.data.frame(peakAnno@anno)
annot_promo <- subset(annot, annotation == "Promoter")
#annot_promo_genes <- annot_promo %>% distinct(ENSEMBL, .keep_all = TRUE)

#print(paste("rows in annotated CGIs -", nrow(annot_promo)))
#print(head(annot_promo))

plotCGI <- merge(annot_promo, CGIdf)
plotCGI <- subset(plotCGI, select = c(seqnames, start, end, ENSEMBL, SYMBOL, pvalue, qvalue, meth.diff))
plotCGI <- merge(KDM5C_binding, plotCGI)
print(paste("rows in CGI merged -", nrow(plotCGI)))
print(head(plotCGI))

write.table(plotCGI, snakemake@output[[2]], sep = ",",  row.names = FALSE)


### volcano plot
#cut off percent methylation difference (10)
l2fcco <- 10
#n = position in snakefile for output volcano plot
#title = name of the dataset for the graph
#restab = methylation results table for plotting
tissue_volcano <- function(n, title, restab){
	#restab <- methyl
	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$keyvals <- ifelse(restab$meth.diff >= l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Bound", 'deeppink',
						ifelse(restab$meth.diff >= l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Unbound", 'palevioletred1',
                        	ifelse(restab$meth.diff <= -l2fcco & restab$qvalue <= 0.01, 'aquamarine',
								ifelse(restab$meth.diff <= -l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Bound", 'deeppink',
									"gray38"))))

	#print(head(restab))

	restab$names <- ifelse(restab$keyvals == 'deeppink', 'KDM5C-bound', 
                         ifelse(restab$keyvals == 'gray38', 'non-significant',
						 	ifelse(restab$keyvals == 'aquamarine', 'hypomethylated',
								ifelse(restab$keyvals == 'palevioletred1', 'hypermethylated',
                            		"uh oh"))))



	#for the figure legend
	keyvals <- restab$keyvals
	names(keyvals) <- restab$names 

	# #Make the axis range narrower so that there is space to see the data
	# neglog10 <- function(x){
	#   a = log10(x)
	#   return(-a)
	# }

	# l2fcmaxie <- 3
	# keyvals.shape <- ifelse(
  	# 	restab$log2FoldChange>l2fcmaxie | neglog10(restab$padj)>10, 9, 16)
	# #if the log2fc or the padj is greater than the maximum axis cut off, set it to the cut off and make the point an outlier shape
	# restab$log2FoldChange[restab$log2FoldChange>l2fcmaxie] <- l2fcmaxie
	# restab$padj[restab$padj< 1e-10] <- 1e-10

	# keyvals.shape[is.na(keyvals.shape)] <- 16
	# names(keyvals.shape)[keyvals.shape == 9] <- 'Outlier'
	# names(keyvals.shape)[keyvals.shape == 16] <- 'Regular'

	#select label
	labels <- c("D1Pas1", "Stra8", "Zar1", "Dazl")

	library(EnhancedVolcano)
	p <- EnhancedVolcano(restab, lab = restab$SYMBOL,
    	                 x = 'meth.diff', y = 'qvalue',
            	         selectLab = labels,
                	     xlim=c(-50,100),
    	                 xlab ="% methylation difference",
        	             ylab = bquote(~-Log[10]~adjusted~italic(P)),
            	         pCutoff = 0.01, FCcutoff = l2fcco,
    	                 title = paste0("WGBS - WT ESC vs exEpiLC ", title),
        	             subtitle = " ",
            	         labSize = 4.0,
                	     colAlpha = .5,
						 pointSize = 3,
                    	 colCustom = keyvals,
						 #shapeCustom = keyvals.shape,
            	         legendPosition = 'right',
                	     gridlines.major = FALSE,
                    	 gridlines.minor = FALSE,
						 drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = 'black',
                	     legendLabSize = 10, legendIconSize = 3.0)

	ggsave(snakemake@output[[n]], plot = p, width = 7, height = 6)


}

tissue_volcano(3, "germline TSS", plotTSS)
tissue_volcano(4, "germline CGI", plotCGI)

