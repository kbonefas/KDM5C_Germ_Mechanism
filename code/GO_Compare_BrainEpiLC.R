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
samples <- c("nESC", "EpiLC", "exEpiLC", "AMY", "HIP")
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
print(head(germDEGs))

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

ggsave(snakemake@output[[9]], p2, width = 6, height = 5.5)





    #ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
##3) Calculate the overlap between groups

#make an upset plot for the overlap
library("UpSetR")
BiocManager::install("UpSetR")

modifiedupset <- function(samplelist){
	upset(fromList(samplelist), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "# of Overlapping Germline DEGs", empty.intersections = "on")
}



#just male samples
pdf(file = snakemake@output[[8]], width = 5.5, height = 5)

upset(fromList(germDEGs), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "Overlapping Germline DEGs", text.scale = 2)

dev.off()


