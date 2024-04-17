#moved calculation of tissue specific genes for brain here to make the code more robust
source("code/utilities/tissueSPECIFICgenes.R") #read in the functions


### run for the samples, make sure order matches snakefile
titles <- c("Amygdala", "Hippocampus")
fisherdfs <- list()
# get the testis-specific DEGs for each samples for later plotting in wt and germline-depleted testis
testisDEGs <- data.frame()
for(i in 1:length(titles)){
	tissDEG <- tissue_genes(i)
	#testis DEGs
	TEdeg <- subset(tissDEG, tissue == "Testis")
	TEdeg <- subset(TEdeg, select = c("ENSEMBL", "SYMBOL"))
	testisDEGs <- rbind(testisDEGs, TEdeg)

	#n = position in snakefile for input DEG table and output volcano plot
	#title = name of the dataset for the graph
	#td = tissue-specific DEG 
	#offset - number of samples to shift the input and output numbers to align correctly
	#bary - y-axis of bar graph
	tissue_plots(i, titles[i], tissDEG, length(titles), 50)


	#calculate the enrichment
	fisherdfs[[i]] <- tissue_fisher(tissDEG, i, titles[i])

}


#save the calculations of tissue number and enrichment

#merge together the tissue count with the fisher results
fisherdfs[[3]] <- countingdf
library(tidyverse)
fisherdfs <- fisherdfs %>% reduce(full_join, by = "tissue")
print(fisherdfs)

write.table(fisherdfs, snakemake@output[[7]], sep = ",", row.names = FALSE)

#saving the testis DEGs
testisDEGs <- unique(testisDEGs)
print(head(testisDEGs))
write.table(testisDEGs, snakemake@output[[8]], sep = ",", row.names = FALSE)



### make a histogram of the total number of tissue specific genes
p <- ggbarplot(countingdf, "tissue", "count", fill = "tissue", color = "tissue", palette = tissue_palette, 
		label = TRUE, lab.pos = "out", lab.col = "black") + rotate_x_text(45) 


pl <- ggpar(p, legend = "none", ylab = "# of genes", xlab = FALSE, title = paste0("Total Tissue Gene #") )
ggsave(snakemake@output[[9]], plot = pl, width = 5, height = 5)

### save tissue-specific DEGs in a dataframe
for(i in 1:length(titles)){
	tissDEG <- tissue_genes(i)
	write.table(tissDEG, snakemake@output[[9+i]], sep = ",", row.names = FALSE)
	}