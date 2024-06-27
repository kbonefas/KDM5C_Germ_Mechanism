## 2024.06.27 Get the transcript start and end site of each gene in the clusters for bedtools plotting of KDM5C binding
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(biomaRt)

genes <- genes(txdb)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#get the list of germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",") 

#attributes = listAttributes(mart, page = "structure")

#get list of genes for each cluster then get the TSS and TES
genelist <- germ$SYMBOL
tss <- getBM(attributes = c("chromosome_name", "transcription_start_site", "transcript_end",
                            "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"), filters = "external_gene_name", values = genelist, mart = mart)

bed <- data.frame(chr = paste0("chr", tss$chromosome_name), start = tss$transcription_start_site - 500, end =  tss$transcription_start_site + 500)

print(head(bed))
write.table(bed, snakemake@output[[1]], sep = "\t", row.names = FALSE, col.names=FALSE, quote = FALSE)
	
