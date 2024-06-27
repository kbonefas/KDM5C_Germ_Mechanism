## 2024.06.27 Get the transcript start and end site of each gene in the clusters for bedtools plotting of KDM5C binding
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(biomaRt)

genes <- genes(txdb)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#get the list of germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",") 

#genelist - list of gene symbols you want the window for
#window - # of bp upstream and downstream
#out - output location in snakemake for bedfile
TSSwindowBED <- function(genelist, window, out){
	#attributes = listAttributes(mart, page = "structure")
	#get the TSS for list of genes
	tss <- getBM(attributes = c("chromosome_name", "transcription_start_site", "transcript_end",
                            "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"), filters = "external_gene_name", values = genelist, mart = mart)

	#add and subtract bp window and format in bed style
	bed <- data.frame(chr = paste0("chr", tss$chromosome_name), start = tss$transcription_start_site - window, end =  tss$transcription_start_site + window)

	write.table(bed, snakemake@output[[out]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)
}

	
WIND <- 1000
TSSwindowBED(germ$SYMBOL, WIND, 1)
testgenes <- c("Dazl", "Cyct", "D1Pas1", "Ddx4")

TSSwindowBED(testgenes, WIND, 2)


#just the brain or EpiLC DEGs
samples <- c("EpiLC", "Amygdala", "Hippocampus")
allDEGs <- c()
for (i in 1:length(samples)){
	germDEGs <- read.csv(snakemake@input[[i+1]], sep = ",")
	allDEGs <- c(allDEGs, germDEGs$SYMBOL)
}

allDEGs <- unique(allDEGs)
print(allDEGs[1:10])
print(length(allDEGs))

TSSwindowBED(allDEGs, WIND, 3)



