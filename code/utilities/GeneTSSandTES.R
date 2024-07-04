#script for getting bed files of the TSS and TES for genes of interest - input is the gtf file used in RNA/ChIPseq quantification

library("rtracklayer")
gtf <- rtracklayer::import('../../gencode.vM21.annotation.gtf') #gtf file used in ChIPseq quantification
gtf_df <- as.data.frame(gtf)
# print(head(gtf_df))

library(dplyr)
library(tidyr)
print(head(gtf_df))
gtf_df <- subset(gtf_df, type == "gene")

#need to switch the start and end if the gene is on the minus strand instead of +
gtf_df$TSS <- ifelse(gtf_df$strand == "+", gtf_df$start, gtf_df$end)
gtf_df$TES <- ifelse(gtf_df$strand == "+", gtf_df$end, gtf_df$start)

#removing the ensembl gene IDs
gtf_df <- pivot_longer(gtf_df, cols =gene_id) %>%
 	separate(value, into = c('ENSEMBL', 'version'), sep = "\\.") %>%
  	select(-name)

	# gtf_df[c('ENSEMBL', 'version')] <- str_split(gtf_df$gene_id, '.')
print('formatted gtf')
print(head(gtf_df))

#get TSS and TES for genes of interest
#goi - genes of interest - vector of ENSEMBL gene IDs 
#n - location of output table in snakefile
geneTSSandTES <- function(goi, n){
	#split the gene_id column into ensembl and variant names

	#regions of interest
	roi <- subset(gtf_df, ENSEMBL %in% goi, select = c(seqnames, TSS, TES))
	print('regions of interest')
	print(tail(roi))

	write.table(roi, snakemake@output[[n]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

}

#goi - genes of interest - vector of ENSEMBL gene IDs 
#window - window of bp around the tss
#n - location of output table in snakefile
geneTSSwindow <- function(goi, n, window){
	#regions of interest
	roi <- subset(gtf_df, ENSEMBL %in% goi, select = c(seqnames, TSS))
	
	#Add and subtract the window from the TSS
	roi$TSS_up <- roi$TSS - window
	roi$TSS_down <- roi$TSS + window

	write.table(roi[,c("seqnames","TSS_up", "TSS_down")], snakemake@output[[n]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

}
