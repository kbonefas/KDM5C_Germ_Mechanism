
### 24.07.26 compare binding locations of KDM5C ChIp-seq peaks and number of germline genes bound for CGI and CGI-free germline genes


# Visualize genomic locations using ChIPseeker
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ggplot2)
library(dplyr)


#get the parameter used for the promoter range
source("code/utilities/parameters.R")


#read in the bedfile

#using WT no KO for both EpiLCs and PNCs
samples <- c("EpiLC", "PNC") #make sure order matches snakefile
#a list to save the plots in
annoplot <- list()
#a list of the peak gene names
peakENSEMBL <- list()

#read in the bed file of CGIs
peak <- readPeakFile(snakemake@input[[1]])
print(head(peak))

#annotate where the peaks are located and plot
peakAnno <- annotatePeak(peak, tssRegion=c(-promorange, promorange), TxDb=txdb, annoDb="org.Mm.eg.db")

#subset for just the promoter peaks
annot <- as.data.frame(peakAnno@anno)
promo <- subset(annot, annotation == "Promoter")

print("promo")
print(head(promo))
	
promo_genes <- promo %>% distinct(ENSEMBL, .keep_all = TRUE) #keep only distinct annotations (in case 2 CGIs are within 1000bp)
CGI_ENSEMBL <- unique(promo_genes$ENSEMBL)
print("promo_genes")
print(head(promo_genes))

print("CGI_ENSEMBL")
print(CGI_ENSEMBL[1:5])


#read in the germline genes list with KDM5C binding at promoter
germ <- read.csv(snakemake@input[[2]], sep = ",")

### save a bed file of the CGI regions that are at germline gene promoters
germ_promo_CGI_bed <- promo[promo$ENSEMBL %in% germ$ENSEMBL,]

write.table(germ_promo_CGI_bed[1:3], snakemake@output[[1]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


germ$Promo_CGI <- ifelse(germ$ENSEMBL %in% CGI_ENSEMBL, "CGI", "no")
write.table(germ, snakemake@output[[2]], sep = ",",  row.names = FALSE)


#get the percentage of all germline genes with promoter CGIs and those bound or unbound by KDM5C
plotdf <- data.frame(Kdm5c_binding = c(rep("KDM5C-Bound", 2), rep("KDM5C-Unbound", 2), rep("All germ", 2)),  CpG_island = rep(c("no", "CGI"), 3))
plotdf
# Kdm5c-Bound	no
# Kdm5c-Bound	CGI
# Kdm5c-Unbound	no
# Kdm5c-Unbound	CGI
# All germ		no
# All germ		CGI

plotdf$CGI_Count <- c(nrow(subset(germ, germ$KDM5C_binding == "Bound" & germ$Promo_CGI == "no")), 
						nrow(subset(germ, germ$KDM5C_binding == "Bound" & germ$Promo_CGI == "CGI")),
						nrow(subset(germ, germ$KDM5C_binding == "Unbound" & germ$Promo_CGI == "no")),
						nrow(subset(germ, germ$KDM5C_binding == "Unbound" & germ$Promo_CGI == "CGI")),
						nrow(subset(germ, germ$Promo_CGI == "no")),
						nrow(subset(germ, germ$Promo_CGI == "CGI")))


#percentage of each category
plotdf$CGIPercent_plot <- ifelse(plotdf$Kdm5c_binding == "KDM5C-Bound", plotdf$CGI_Count/sum(subset(plotdf, plotdf$Kdm5c_binding == "KDM5C-Bound")$CGI_Count) * 100, 
								ifelse(plotdf$Kdm5c_binding == "KDM5C-Unbound", plotdf$CGI_Count/sum(subset(plotdf, plotdf$Kdm5c_binding == "KDM5C-Unbound")$CGI_Count) * 100,
								ifelse(plotdf$Kdm5c_binding == "All germ", plotdf$CGI_Count/nrow(germ) * 100, 0)))

plotdf$Percent <- as.integer(round(plotdf$CGIPercent_plot))
plotdf


#plot the results in a bar graph
#set plotting order
plotdf$Kdm5c_binding <- factor(plotdf$Kdm5c_binding, levels = c("All germ", "KDM5C-Bound", "KDM5C-Unbound"))

library("ggpubr")
q <- ggbarplot(plotdf, "Kdm5c_binding", "Percent", fill = "CpG_island", color = "CpG_island", palette = c("no" = "aquamarine4", "CGI" = "aquamarine3"),
		title = "CpG islands at germline gene promoters", label = TRUE, lab.col = "black", lab.vjust = 1, xlab = " ", ylab = "% of genes", orientation = "vert") 

ggsave(snakemake@output[[3]], q, width = 5, height = 4)


#gene ontology of CGI germline genes vs non
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# make a list with all the germline genes
#sample names, make sure order matches snakefile input
samples <- c("CGI", "no")
germDEGs <- list()

for(i in 1:length(samples)){
	#get the ENSEMBL names of CGI vs non CGI promoters
	germDEGs[[i]] <- germ[germ$Promo_CGI == samples[i], ][,1]
}

names(germDEGs) <- samples


#now run the gene ontology comparison
ck <- compareCluster(geneCluster = germDEGs, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck) 

write.table(ck, snakemake@output[[4]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[5]], plot = dotplot(ck, showCategory = 10), width = 5, height = 9)



### Save the gene names for HOMER motif analysis
write.table(subset(germ, Promo_CGI == "CGI")[,"SYMBOL"], snakemake@output[[6]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(subset(germ, Promo_CGI == "no")[,"SYMBOL"], snakemake@output[[7]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
