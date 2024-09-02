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

#function to annotate TSS df
anno_TSS <- function(TSSdf){
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
	return(plotTSS)
}


## TSS WT ESC vs EpiLCs
ESCvsEpiLC_TSS <- anno_TSS(read.csv(snakemake@input[[2]], sep = ","))
write.table(ESCvsEpiLC_TSS, snakemake@output[[1]], sep = ",", row.names = FALSE)

## TSS WT vs KO
WTvsKO_TSS <- anno_TSS(read.csv(snakemake@input[[4]], sep = ","))
write.table(WTvsKO_TSS, snakemake@output[[5]], sep = ",", row.names = FALSE)


### TSS volcano plot
#cut off percent methylation difference 
l2fcco <- 25

### highlighting CGI vs not genes
# n is the location of the output file
TSS_volcano <- function(n, TITLE, restab, maxy){

	#restab <- methyl
	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$keyvals <- ifelse(restab$qvalue > 0.01 & restab$Promo_CGI == "CGI", 'black',
						ifelse(restab$Promo_CGI == "CGI" & restab$meth.diff < l2fcco & restab$meth.diff > -l2fcco, 'black',
							ifelse(restab$qvalue <= 0.01 & restab$Promo_CGI == "CGI", 'red',
								ifelse(restab$meth.diff >= l2fcco & restab$qvalue <= 0.01 & restab$Promo_CGI == "no", 'palevioletred1',
                        			ifelse(restab$meth.diff <= -l2fcco & restab$qvalue <= 0.01 & restab$Promo_CGI == "no", 'lightseagreen',
										"gray38")))))

	#print(head(restab))

	restab$names <- ifelse(restab$keyvals == 'red', 'promoter CGI', 
                         ifelse(restab$keyvals == 'gray38', 'n.s. non-CGI',
						 	ifelse(restab$keyvals == 'black', 'n.s. promoter CGI',
								ifelse(restab$keyvals == 'lightseagreen', 'non-CGI hypomethylated',
									ifelse(restab$keyvals == 'palevioletred1', 'non-CGI hypermethylated',
                            			"uh oh")))))



	#for the figure legend
	keyvals <- restab$keyvals
	names(keyvals) <- restab$names 

	#Make the axis range narrower so that there is space to see the data
	neglog10 <- function(x){
	  a = log10(x)
	  return(-a)
	}


	
	keyvals.shape <- ifelse(neglog10(restab$qvalue) > maxy, 9, 16)
	#if the log2fc or the padj is greater than the maximum axis cut off, set it to the cut off and make the point an outlier shape
	restab$qvalue[restab$qvalue< 10^(-maxy)] <- 10^(-maxy)

	keyvals.shape[is.na(keyvals.shape)] <- 16
	names(keyvals.shape)[keyvals.shape == 9] <- 'Outlier'
	names(keyvals.shape)[keyvals.shape == 16] <- 'Regular'

	
	#labdesc <- ifelse(labels == "yes", restab$SYMBOL, ifelse(labels == "no", NA, "oop"))
	labels <- c("D1Pas1", "Cyct", "Ddx4", "Stra8", "Zar1", "Dazl", "Naa11", "Tex15", "Tex14", "Tex11", "Rpl10l")

	library(EnhancedVolcano)
	p <- EnhancedVolcano(restab, lab = restab$SYMBOL,
    	                 x = 'meth.diff', y = 'qvalue',
            	         selectLab = labels,
                	     xlim=c(-100,100),
    	                 xlab ="% methylation difference",
        	             ylab = bquote(~-Log[10]~qvalue),
            	         pCutoff = 0.01, FCcutoff = l2fcco,
    	                 title = TITLE,
        	             subtitle = " ",
            	         labSize = 5.5,
                	     colAlpha = .6,
						 pointSize = 3,
                    	 colCustom = keyvals,
						 shapeCustom = keyvals.shape,
            	         legendPosition = 'right',
                	     gridlines.major = FALSE,
                    	 gridlines.minor = FALSE,
						 drawConnectors = TRUE, widthConnectors = .5, colConnectors = 'black',
                	     legendLabSize = 12, legendIconSize = 3.0)

	ggsave(snakemake@output[[n]], plot = p, width = 9, height = 7)


}

###for kdm5c binding
TSS_volcano_5C <- function(n, TITLE, restab, maxy){

	#restab <- methyl
	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$keyvals <- ifelse(restab$qvalue > 0.01 & restab$KDM5C_binding == "Bound", 'black',
						ifelse(restab$KDM5C_binding == "Bound" & restab$meth.diff < l2fcco & restab$meth.diff > -l2fcco, 'black',
							ifelse(restab$qvalue <= 0.01 & restab$KDM5C_binding == "Bound", 'blue',
								ifelse(restab$meth.diff >= l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Unbound", 'palevioletred1',
                        			ifelse(restab$meth.diff <= -l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Unbound", 'lightseagreen',
										"gray38")))))

	#print(head(restab))

	restab$names <- ifelse(restab$keyvals == 'blue', 'KDM5C-bound', 
                         ifelse(restab$keyvals == 'gray38', 'n.s. KDM5C-unbound',
						 	ifelse(restab$keyvals == 'black', 'n.s. KDM5C-bound',
								ifelse(restab$keyvals == 'lightseagreen', 'KDM5C-unbound',
									ifelse(restab$keyvals == 'palevioletred1', 'KDM5C-unbound',
                            			"uh oh")))))



	#for the figure legend
	keyvals <- restab$keyvals
	names(keyvals) <- restab$names 

	#Make the axis range narrower so that there is space to see the data
	neglog10 <- function(x){
	  a = log10(x)
	  return(-a)
	}


	
	keyvals.shape <- ifelse(neglog10(restab$qvalue) > maxy, 9, 16)
	#if the log2fc or the padj is greater than the maximum axis cut off, set it to the cut off and make the point an outlier shape
	restab$qvalue[restab$qvalue< 10^(-maxy)] <- 10^(-maxy)

	keyvals.shape[is.na(keyvals.shape)] <- 16
	names(keyvals.shape)[keyvals.shape == 9] <- 'Outlier'
	names(keyvals.shape)[keyvals.shape == 16] <- 'Regular'

	
	#labdesc <- ifelse(labels == "yes", restab$SYMBOL, ifelse(labels == "no", NA, "oop"))

	library(EnhancedVolcano)
	p <- EnhancedVolcano(restab, lab = restab$SYMBOL,
    	                 x = 'meth.diff', y = 'qvalue',
            	         #selectLab = labels,
                	     xlim=c(-100,100),
    	                 xlab ="% methylation difference",
        	             ylab = bquote(~-Log[10]~qvalue),
            	         pCutoff = 0.01, FCcutoff = l2fcco,
    	                 title = TITLE,
        	             subtitle = " ",
            	         labSize = 4.0,
                	     colAlpha = .6,
						 pointSize = 3,
                    	 colCustom = keyvals,
						 shapeCustom = keyvals.shape,
            	         legendPosition = 'right',
                	     gridlines.major = FALSE,
                    	 gridlines.minor = FALSE,
						 drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = 'black',
                	     legendLabSize = 10, legendIconSize = 3.0)

	ggsave(snakemake@output[[n]], plot = p, width = 9, height = 7)


}

TSS_volcano(3, "ESC vs exEpiLC - germline TSS", ESCvsEpiLC_TSS, 200)
TSS_volcano(7, "WT vs KO exEpiLC - germline TSS", WTvsKO_TSS, 50)
TSS_volcano_5C(9, "WT vs KO exEpiLC - germline TSS", WTvsKO_TSS, 50)




## germline CGIs
anno_CGI <- function(CGIdf){
	
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

	return(plotCGI)

}

ESCvsEpiLC_CGI <- anno_CGI(read.csv(snakemake@input[[3]], sep = ","))
write.table(ESCvsEpiLC_CGI, snakemake@output[[2]], sep = ",",  row.names = FALSE)


WTvsKO_CGI <- anno_CGI(read.csv(snakemake@input[[5]], sep = ","))
write.table(WTvsKO_CGI, snakemake@output[[6]], sep = ",",  row.names = FALSE)



### Highlighting KDM5C-bound vs unbound CGIs
#n = position in snakefile for output volcano plot
#title = name of the dataset for the graph
#restab = methylation results table for plotting
CGI_volcano <- function(n, TITLE, restab){
	#restab <- methyl
	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$keyvals <- ifelse(restab$qvalue <= 0.01 & restab$KDM5C_binding == "Bound"  & restab$meth.diff >= l2fcco, 'mediumslateblue',
						ifelse(restab$qvalue > 0.01 & restab$KDM5C_binding == "Bound" & restab$meth.diff < l2fcco, 'black',
							ifelse(restab$meth.diff >= l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Unbound", 'palevioletred1',
                        		ifelse(restab$meth.diff <= -l2fcco & restab$qvalue <= 0.01 & restab$KDM5C_binding == "Unbound", 'lightseagreen',
									"gray38"))))

	#print(head(restab))

	restab$names <- ifelse(restab$keyvals == 'mediumslateblue', 'KDM5C-bound', 
                         ifelse(restab$keyvals == 'gray38', 'n.s.',
						 	ifelse(restab$keyvals == 'black', 'n.s. KDM5C-bound',
								ifelse(restab$keyvals == 'lightseagreen', 'hypomethylated',
									ifelse(restab$keyvals == 'palevioletred1', 'hypermethylated',
                            			"uh oh")))))



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
	#labels <- c("D1Pas1", "Cyct", "Stra8", "Zar1", "Dazl", "Naa11", "Tex15", "Tex14", "Tex11", "Rpl10l")

	library(EnhancedVolcano)
	p <- EnhancedVolcano(restab, lab = restab$SYMBOL,
    	                 x = 'meth.diff', y = 'qvalue',
            	         #selectLab = labels,
                	     xlim=c(-100,100),
    	                 xlab ="% methylation difference",
        	             ylab = bquote(~-Log[10]~adjusted~italic(P)),
            	         pCutoff = 0.01, FCcutoff = l2fcco,
    	                 title = TITLE,
        	             subtitle = " ",
            	         labSize = 4.0,
                	     colAlpha = .6,
						 pointSize = 3,
                    	 colCustom = keyvals,
						 #shapeCustom = keyvals.shape,
            	         legendPosition = 'right',
                	     gridlines.major = FALSE,
                    	 gridlines.minor = FALSE,
						 drawConnectors = FALSE, widthConnectors = 1.0, colConnectors = 'black',
                	     legendLabSize = 10, legendIconSize = 3.0)

	ggsave(snakemake@output[[n]], plot = p, width = 7, height = 5)


}

CGI_volcano(4, "ESC vs exEpiLC - germline CGI", ESCvsEpiLC_CGI)
CGI_volcano(8, "WT vs KO - germline CGI", WTvsKO_CGI)


# #write bedfile with the coordinates
# getbed <- function(df, kdm5c, n){
# 	bed <- subset(df, KDM5C_binding == kdm5c, select = c("seqnames", "start", "end"))	
# 	write.table(bed, snakemake@output[[n]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)
# }

# getbed(plotTSS, "Bound", 5)
# getbed(plotTSS, "Unbound", 6)
# getbed(plotCGI, "Bound", 7)
# getbed(plotCGI, "Unbound", 8)


