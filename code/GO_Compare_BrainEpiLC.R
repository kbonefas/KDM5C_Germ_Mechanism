#24.04.09 Gene ontology comparison between male Brain DEGs and EpiLC germline genes to demonstrate the different types of genes dysregulated
#updated 25.09.01 to include all RNAseq datasets
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# make a list with all the germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",", header = TRUE)
print("germ")
print(head(germ))

#sample names, make sure order matches snakefile input
samples <- c("nESC", "EpiLC", "exEpiLC", "NPC", "AMY", "HIP")
germDEGs <- list()

for(i in 1:length(samples)){
	DEGs <- read.csv(snakemake@input[[i+1]], sep = ",", header = TRUE)
	print(head(DEGs))

	#subset for germline DEGs
	germDEGs1 <- subset(DEGs, DEGs$ENSEMBL %in% germ$ENSEMBL)
	
	#save the germ DEGs
	write.table(germDEGs1, snakemake@output[[i]], sep = ",", row.names = FALSE)

	#get the ensembl names, put in the position of the list
	germDEGs[[i]] <- germDEGs1[,1]
}

names(germDEGs) <- samples


#now run the gene ontology comparison
ck <- compareCluster(geneCluster = germDEGs, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck)

write.table(ck, snakemake@output[[length(samples) + 1]], row.names = FALSE, sep = ",")

p <- dotplot(ck, size = "Count") +
  theme(axis.text.y = element_text(size=8)) +
  scale_color_gradient(low = "blue3", high = "red")

ggsave(snakemake@output[[length(samples) + 2]], p, width = 6, height = 5.5)


### simplify the clusters 
simple <- simplify(ck)
p2 <- dotplot(simple, size = "Count") +
  theme(axis.text.y = element_text(size=8)) +
  scale_color_gradient(low = "blue3", high = "red")

ggsave(snakemake@output[[length(samples) + 4]], p2, width = 6, height = 5.5)





    #ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
##3) Calculate the overlap between groups

#make an upset plot for the overlap
library("UpSetR")
# BiocManager::install("UpSetR")

modifiedupset <- function(samplelist){
	upset(fromList(samplelist), sets.x.label = "# Germline DEGs", mainbar.y.label = "# of Overlapping Germline DEGs")
}


#print(germ in)
#change the plotting order
differentiation <- c( "AMY", "HIP", "NPC", "exEpiLC", "EpiLC", "nESC")
germDEGs2 <- germDEGs[differentiation]
print(head(germDEGs2))

pdf(file = snakemake@output[[length(samples) + 3]], width = 8, height = 5.5)

upset(fromList(germDEGs2), order.by = "freq",  sets.x.label = "# germline DEGs", mainbar.y.label = "# in group", text.scale = 2, sets = differentiation, mb.ratio = c(0.7, 0.3), keep.order = TRUE)

dev.off()


# 2026.07.07 percentage of all germline genes that are DEGs at any point
allgermDEGs<- unique(unlist(germDEGs))
print(paste0("Number of germline genes that are DEGs at any point: ", length(allgermDEGs)))

#subset germ genes that are not DEGs
nonDEGs <- subset(germ, !(germ$ENSEMBL %in% allgermDEGs))
print(paste0("Number of germline genes that are not DEGs at any point: ", nrow(nonDEGs)))

print(paste0("Number of germline genes: ", nrow(germ)))


library(ggplot2)

# Create Data
data <- data.frame(
  group=c("DEG", "non-DEG"),
  value=c(length(allgermDEGs)/nrow(germ), nrow(nonDEGs)/nrow(germ))
)
print(data)
data$percent <- round(data$value * 100, 2)
data$labels <- paste0(data$percent, "%")
# Basic piechart
source("code/utilities/colorpalettes.R")
p <- ggplot(data, aes(x="", y=percent, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  
  theme_void() + # remove background, grid, numeric labels

  geom_text(aes(label = labels), position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_manual(values=c("DEG" = germcolor, "non-DEG" = "#194219"))

pdf(file = snakemake@output[[length(samples) + 5]], width = 5.5, height = 5)

p

dev.off()