
### 24.04.26 compare binding locations of KDM5C ChIp-seq peaks and number of germline genes bound


# Visualize genomic locations using ChIPseeker
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ggplot2)
library(dplyr)



#get the parameter used for the promoter range
source("code/utilities/parameters.R")

#using WT no KO for both EpiLCs and PNCs
samples <- c("EpiLC", "PNC") #make sure order matches snakefile
#a list to save the plots in
annoplot <- list()
#a list of the peak gene names
peakENSEMBL <- list()


for (i in 1:length(samples)){
	#make a GRanges object from your .bed file of peak locations
	peak <- readPeakFile(snakemake@input[[i]])
	print(head(peak))

	#annotate where the peaks are located and plot
	
	peakAnno <- annotatePeak(peak, tssRegion=c(-promorange, promorange),
                         TxDb=txdb, annoDb="org.Mm.eg.db")


	#subset for just the promoter peaks
	annot_KDM5C <- as.data.frame(peakAnno@anno)
	annot_KDM5C.promo <- subset(annot_KDM5C, annotation == "Promoter")
	KDM5C.promo_genes <- annot_KDM5C.promo %>% distinct(ENSEMBL, .keep_all = TRUE) #keep only distinct annotations (in case 2 peaks are within 1000bp)
	#print(head(KDM5C.promo_genes))
	write.table(KDM5C.promo_genes, snakemake@output[[i]], sep=",", row.names = FALSE)



	peakENSEMBL[[i]] <- unique(KDM5C.promo_genes$ENSEMBL)

	annoplot[[i]] <- plotAnnoBar(peakAnno, title = paste(samples[i], "-", format(nrow(annot_KDM5C), big.mark = ","), "peaks"))


}

library("ggplot2")
library("gridExtra")
library("ggpubr")
ggsave(snakemake@output[[3]], grid.arrange(grobs = annoplot, nrow = 2), width = 6, height = 4)



#See which genes that kdm5c is bound to in epilcs and pncs overlap with DEGs
#put the name of the peak samples with their ensembls
names(peakENSEMBL) <- samples

#read in the germline DEGs
DEGsample <- c("Amygdala", "Hippocampus", "EpiLC")
germDEGs <- list()

for(i in 1:length(DEGsample)){
	DEGs <- read.csv(snakemake@input[[i+length(samples)]], sep = ",", header = TRUE)
	#get the ensembl names, put in the position of the list
	germDEGs[[i]] <- DEGs[,1]
}

names(germDEGs) <- DEGsample
print(head(germDEGs))

#make a new category that is all of the brain DEGs together
germDEGs[["Brain"]] <- unique(c(germDEGs[["Amygdala"]], germDEGs[["Hippocampus"]]))

#get the genes that are in both EpiLcs and brain and unique
sharedgerm <- intersect(germDEGs[["Brain"]], germDEGs[["EpiLC"]])
EpiLCgerm <- setdiff(germDEGs[["EpiLC"]], germDEGs[["Brain"]])
Braingerm <- setdiff(germDEGs[["Brain"]], germDEGs[["EpiLC"]])

germCategory <- c("Common", "EpiLC Only", "Brain Only")
germGroups <- list(sharedgerm, EpiLCgerm, Braingerm)
names(germGroups) <- germCategory

#count the number of promoter peaks genes that fall into each category
	#dataframe for each peakgroup (EpiLC or PNC)
	#number kdm5c bound
	#number kdm5c unbound
qlotdf <- data.frame() #dataframe of counts for plotting

genestatus <- data.frame() #dataframe of which gene is which for motif analysis

#for EpiLc and PNC
for (k in samples){
	for (s in germCategory){
	#read in the dataframe
	a <- germGroups[[s]]
	KDM5C_bound <- a[a %in% peakENSEMBL[[k]]]
	KDM5C_unbound <- a[!a %in% peakENSEMBL[[k]]]
	df <- data.frame(ChIP =  rep(k, 2), GeneCategory = rep(s, 2), KDM5C_binding = c("Bound", "Unbound"), Count = c(length(KDM5C_bound), length(KDM5C_unbound)))
	qlotdf <- rbind(qlotdf, df)

	if(k == "EpiLC"){
		genestatus_df <- data.frame(ENSEMBL = c(KDM5C_bound, KDM5C_unbound), DEG_Type = rep(s, length(KDM5C_bound)+length(KDM5C_unbound)) , KDM5C_binding = c(rep("Bound", length(KDM5C_bound)), rep("Unbound", length(KDM5C_unbound))))
		genestatus <- rbind(genestatus_df, genestatus)
	} else {
		next
	}


	}

}


#plotting data frame
print(qlotdf)
#save the plots in an empty list
plots <- list()

for (k in 1:length(samples)){
	plotdf <- subset(qlotdf, ChIP == samples[k])
	#set the plotting order
	plotdf$GeneCategory <- factor(plotdf$GeneCategory, levels = c("EpiLC Only", "Common", "Brain Only"))

	library("ggpubr")
	q <- ggbarplot(plotdf, "GeneCategory", "Count",
  		fill = "KDM5C_binding", color = "KDM5C_binding", palette = c("Bound" = "darkorange1", "Unbound" = "lightskyblue4" ),
  		title = paste(samples[k], "peaks") , label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "RNA-seq germline DEG", ylab = "# of germline DEGs", orientation = "vert") 

	plots[[k]] <- q

}

ggsave(snakemake@output[[4]], grid.arrange(grobs = plots, ncol = 2), width = 7.5, height = 3.5)


### Which germline genes are bound by KDM5C
germ <- read.csv(snakemake@input[[6]], sep = ",")
# print("germ genes")
# # print(head(germ))

# print("genestatus")
genestatus <- merge(genestatus, germ, by = "ENSEMBL")
# print(head(genestatus))

#save all of the bound and unbound DEGs in EpiLCs
write.table(genestatus, snakemake@output[[5]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#save list of genes for HOMER
write.table(subset(genestatus, KDM5C_binding == "Bound")[,"SYMBOL"], snakemake@output[[6]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(subset(genestatus, KDM5C_binding == "Unbound")[,"SYMBOL"], snakemake@output[[7]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#all germline genes, which ones are kdm5c bound to in EpiLCs
germ$KDM5C_binding <- ifelse(germ$ENSEMBL %in% peakENSEMBL[["EpiLC"]], "Bound", "Unbound")
print("all germ status")
print(head(germ))

#save table
write.table(germ, snakemake@output[[8]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#save list of genes for HOMER
write.table(subset(germ, KDM5C_binding == "Bound")[,"SYMBOL"], snakemake@output[[9]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(subset(germ, KDM5C_binding == "Unbound")[,"SYMBOL"], snakemake@output[[10]], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#make eulerr of kdm5c binding at all germline genes
library("eulerr")
source("code/utilities/colorpalettes.R")

#genes bound in epilc that are not germline genes
allepi <- peakENSEMBL[["EpiLC"]]
epilc_nogerm <- allepi[!allepi %in% germ$ENSEMBL]


together <- euler(c("Bound_EpiLCs" = length(epilc_nogerm), "Germline" = nrow(subset(germ, KDM5C_binding == "Unbound")), "Bound_EpiLCs&Germline" = nrow(subset(germ, KDM5C_binding == "Bound"))))

p <- plot(together, quantities = TRUE, labels = list(font = 4), fills = c(EpiLC_XY_KO, germcolor))

library("ggplot2")
ggsave(snakemake@output[[11]], plot = p, width = 4, height = 3)
