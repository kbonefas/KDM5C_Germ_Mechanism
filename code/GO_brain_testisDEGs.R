#24.04.18 Gene ontology of 5CKO testis DEGs from amygdala and hippocampus
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)


#read in the testis DEGs
testisDEGs <- read.csv(snakemake@input[[1]], sep = ",", header = TRUE)
	#first column is the ENSEMBL IDs

#choose what organism and what type of ontology (biological process)
ego <- enrichGO(testisDEGs[,1], keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
print(ego[, c("ID", "Description", "p.adjust")])

q <- dotplot(ego, showCategory=15)


ggsave(snakemake@output[[1]], plot = q, width = 5.5, height = 6)


#ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
