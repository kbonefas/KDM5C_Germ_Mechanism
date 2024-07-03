#24.04.30 Gene ontology comparison between KDM5C-bound promoters in EpiLCs and PNCs
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

ck <- simplify(ck)
write.table(ck, snakemake@output[[1]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[2]], plot = dotplot(ck, size = "Count"), width = 9, height = 5.5)


#gene ontology for genes only bound in epilc, genes only bound in PNC, and genes bound in both
#genes in all groups
all <- Reduce(intersect, promogenes)
print(head(all))

EpiLC_unique <- promogenes[["EpiLC"]][!promogenes[["EpiLC"]] %in% all]
PNC_unique <- promogenes[["PNC"]][!promogenes[["PNC"]] %in% all]
print(head(EpiLC_unique))

promo_compare <- list(all, EpiLC_unique, PNC_unique)
names(promo_compare) <- c("Shared", "EpiLC only", "PNC only")


#now run the gene ontology comparison
ck2 <- compareCluster(geneCluster = promo_compare, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck2) 
write.table(ck2, snakemake@output[[3]], row.names = FALSE, sep = ",")
ggsave(snakemake@output[[4]], plot = dotplot(ck2, size = "Count"), width = 5.5, height = 5)


#get the overlap between promoters (euler plot)
library("eulerr")
source("code/utilities/colorpalettes.R")
together <- euler(c("EpiLC" = length(EpiLC_unique), "PNC" = length(PNC_unique), "EpiLC&PNC" = length(all)))

p <- plot(together, quantities = TRUE, labels = list(font = 4), fills = c(EpiLC_XY_KO, HIPKO))

library("ggplot2")
ggsave(snakemake@output[[5]], plot = p, width = 3, height = 3)

