#24.04.09 Gene ontology comparison between male Brain DEGs and EpiLC germline genes to demonstrate the different types of genes dysregulated
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# make a list with all the germline genes

#sample names, make sure order matches snakefile input
samples <- c("EpiLC", "PNC")
promogenes <- list()

for(i in 1:length(samples)){
	genes <- read.csv(snakemake@input[[i]], sep = ",", header = TRUE)
	print(head(genes))
	#get the ensembl names, put in the position of the list
	promogenes[[i]] <- genes$ENSEMBL
}

names(promogenes) <- samples


#now run the gene ontology comparison
ck <- compareCluster(geneCluster = promogenes, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck) 

write.table(ck, snakemake@output[[1]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[2]], plot = dotplot(ck, size = "Count"), width = 9, height = 5.5)


    #ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
