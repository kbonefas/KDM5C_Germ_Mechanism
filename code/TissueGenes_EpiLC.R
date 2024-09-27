#24.03.06 - generate graphs for tissue specific genes in EpiLCs and P6 brain

source("code/utilities/tissueSPECIFICgenes.R") #read in the functions


### run for the samples, make sure order matches snakefile
titles <- c("EpiLC XY 5CKO")
fisherdfs <- list()

for(i in 1:length(titles)){
	tissDEG <- tissue_genes(i)

	#n = position in snakefile for input DEG table and output volcano plot
	#title = name of the dataset for the graph
	#td = tissue-specific DEG 
	#offset - number of samples to shift the input and output numbers to align correctly
	#bary - y-axis of bar graph
	bary <- ifelse(i <= 3, 100, 20) #if it's less than or equal to 3 it's the epilc samples which need a higher y-axis
	tissue_plots(i, titles[i], tissDEG, length(titles), bary)

	fisherdfs[[i]] <- tissue_fisher(tissDEG, i, titles[i])

}

#merge together the tissue count with the fisher results
fisherdfs[[length(titles)+1]] <- countingdf
library(tidyverse)
fisherdfs <- fisherdfs %>% reduce(full_join, by = "tissue")
print(fisherdfs)

write.table(fisherdfs, snakemake@output[[4]], sep = ",", row.names = FALSE)


### save tissue-specific DEGs in a dataframe
for(i in 1:length(titles)){
	tissDEG <- tissue_genes(i)
	#all DEGs
	DEGs <- read.csv(snakemake@input[[i]], sep =",")
	# print(paste(n, "tissue DEGs:", nrow(tissue_DEGs)))
	print(paste("all DEGs", nrow(DEGs)))
	print(tail(DEGs))
	#all DEGs together annotated
	annot <- merge(DEGs, tissDEG, by = "ENSEMBL", all.x = T)
	annot <- annot[!duplicated(annot), ]
	print(paste("merge DEGs", nrow(annot)))
	print(tail(annot))
	#annot <- subset(annot, select = c("ENSEMBL", "baseMean.x", "log2FoldChange.x", "lfcSE.x ", "pvalue.x", "padj.x", "Direction.x", "SYMBOL.x", "tissue"))
	write.table(annot, snakemake@output[[i+4]], sep = ",", row.names = FALSE)
	}