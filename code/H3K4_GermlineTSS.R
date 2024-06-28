## 2024.06.27 Get the transcript start and end site of each gene in the clusters for bedtools plotting of KDM5C binding
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(biomaRt)

genes <- genes(txdb)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#get the list of germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",") 

#genelist - list of ensembl gene IDs you want the window for
#window - # of bp upstream and downstream
#out - output location in snakemake for bedfile
TSSwindowBED <- function(genelist, window, out){
	#listAttributes(mart, page = "structure")
	#get the TSS for list of genes
	tss <- getBM(attributes = c("chromosome_name", "transcription_start_site", "transcript_end",
                            "ensembl_gene_id", "start_position", "end_position", "ensembl_transcript_id", "external_gene_name"), filters = "ensembl_gene_id", values = genelist, mart = mart)

	#add and subtract bp window and format in bed style
	print(tss)
	bed <- unique(data.frame(chr = paste0("chr", tss$chromosome_name), start = tss$transcription_start_site - window, end =  tss$transcription_start_site + window))
	print(bed)
	write.table(bed, snakemake@output[[out]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

	bed2 <- unique(data.frame(chr = paste0("chr", tss$chromosome_name), start = tss$start_position, end =  tss$end_position))
	print(bed2)
	write.table(bed2, snakemake@output[[out+1]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

}

	
WIND <- 1000
#TSSwindowBED(germ$ENSEMBL, WIND, 1)
testgenes <- c("ENSMUSG00000010592", "ENSMUSG00000056436", "ENSMUSG00000021758")

TSSwindowBED(testgenes, WIND, 1)


#just the brain or EpiLC DEGs
samples <- c("EpiLC", "Amygdala", "Hippocampus")
allDEGs <- c()
for (i in 1:length(samples)){
	germDEGs <- read.csv(snakemake@input[[i+1]], sep = ",")
	allDEGs <- c(allDEGs, germDEGs$ENSEMBL)
}

allDEGs <- unique(allDEGs)
print(allDEGs[1:10])
print(length(allDEGs))

#TSSwindowBED(allDEGs, WIND, 5)



