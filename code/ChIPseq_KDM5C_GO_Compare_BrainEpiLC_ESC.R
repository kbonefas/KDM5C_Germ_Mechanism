#24.04.30 Gene ontology comparison between KDM5C-bound promoters in EpiLCs and PNCs
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# make a list with all the germline genes

#sample names, make sure order matches snakefile input
samples <- c("pESC", "EpiLC", "PNC")
promogenes <- list()

for(i in 1:length(samples)){
	genes <- read.csv(snakemake@input[[i]], sep = ",", header = TRUE)
	print(head(genes))
	#get the ensembl names, put in the position of the list
	promogenes[[i]] <- genes$ENSEMBL
	print(paste0("number of genes in ", samples[i], " = ", length(promogenes[[i]])))
}

names(promogenes) <- samples


#now run the gene ontology comparison
ck <- compareCluster(geneCluster = promogenes, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck) 

ck <- simplify(ck)
write.table(ck, snakemake@output[[1]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[2]], plot = dotplot(ck, size = "Count"), width = 5, height = 10)


#gene ontology for genes only bound in epilc, genes only bound in PNC, and genes bound in both
#genes in all groups
all <- Reduce(intersect, promogenes)
print(head(all))

pESC_EpiLC <- setdiff(intersect(promogenes[["pESC"]], promogenes[["EpiLC"]]), all)
pESC_PNC <- setdiff(intersect(promogenes[["pESC"]], promogenes[["PNC"]]), all)
EpiLC_PNC <- setdiff(intersect(promogenes[["EpiLC"]], promogenes[["PNC"]]), all)


pESC_unique <- setdiff(promogenes[["pESC"]], c(all, pESC_EpiLC, pESC_PNC))
EpiLC_unique <- setdiff(promogenes[["EpiLC"]], c(all, pESC_EpiLC, EpiLC_PNC))
PNC_unique <- setdiff(promogenes[["PNC"]], c(all, pESC_PNC, EpiLC_PNC))


#check the overlaps
print(paste0("pESC = ", length(pESC_unique) + length(pESC_EpiLC) + length(pESC_PNC) + length(all)))
print(paste0("EpiLC = ", length(EpiLC_unique) + length(pESC_EpiLC) + length(EpiLC_PNC) + length(all)))
print(paste0("PNC = ", length(PNC_unique) + length(pESC_PNC) + length(EpiLC_PNC) + length(all)))


#list with comparisons
promo_compare <- list(all, pESC_unique, EpiLC_unique, PNC_unique, pESC_EpiLC, pESC_PNC, EpiLC_PNC)
names(promo_compare) <- c("Shared", "pESC only", "EpiLC only", "PNC only", "pESC&EpiLC", "pESC&PNC", "EpiLC&PNC")


#now run the gene ontology comparison
ck2 <- compareCluster(geneCluster = promo_compare, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck2) 
write.table(ck2, snakemake@output[[3]], row.names = FALSE, sep = ",")
ggsave(snakemake@output[[4]], plot = dotplot(ck2, size = "Count"), width = 8, height = 10)

#####################
#get the overlap between promoters (euler plot)
library("eulerr")
source("code/utilities/colorpalettes.R")
together <- euler(c("pESC" = length(pESC_unique), "EpiLC" = length(EpiLC_unique), "PNC" = length(PNC_unique), "pESC&EpiLC" = length(pESC_EpiLC), "pESC&PNC" = length(pESC_PNC), "EpiLC&PNC" = length(EpiLC_PNC), "pESC&EpiLC&PNC" = length(all)))

p <- plot(together, quantities = TRUE, labels = list(font = 4), fills = c(ESC_XY_KO, EpiLC_XY_KO, HIPKO))

library("ggplot2")
ggsave(snakemake@output[[5]], plot = p, width = 3.5, height = 3.5)

