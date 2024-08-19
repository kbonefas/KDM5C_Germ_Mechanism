# 2024.07.15 - heatmap of shared XX and XY DEGs

####1) Shared germline genes of interest
XXvsXY <- read.csv(snakemake@input[[1]], sep = ",")

sharedgroups <- c("All", "XY 5CKO and XX 5CHET", "XY 5CKO and XX 5CKO")
XXvsXY <- XXvsXY[XXvsXY$Sample_DEG %in% sharedgroups,]



###2) get the log2fc from WT for each genotype
# clustering with log2fold change from WT
log2fc <- list()
samples <- c("XY_5CKO", "XX_5CHET",  "XX_5CKO")
for (i in 1:length(samples)){
	#read in the results table
	resdf <- read.csv(snakemake@input[[i+1]], row.names = 1)
	#print(head(resdf))
	#get the ensembl and log2fold change
	sampL2FC <- subset(resdf, row.names(resdf) %in% XXvsXY$ENSEMBL, select = "log2FoldChange")
	names(sampL2FC)[names(sampL2FC) == 'log2FoldChange'] <- samples[i] 

	log2fc[[i]] <- sampL2FC

}


library(dplyr)
log2fc <- bind_cols(log2fc)

#make the row names the gene symbol instead of ensembl
log2fc$ENSEMBL <- row.names(log2fc)
log2fc_SYMBOL <- merge(log2fc, XXvsXY, by = "ENSEMBL")
# print(head(log2fc_SYMBOL))

log2fc <- data.frame(row.names = log2fc_SYMBOL$SYMBOL, log2fc_SYMBOL$XY_5CKO, log2fc_SYMBOL$XX_5CHET, log2fc_SYMBOL$XX_5CKO)
colnames(log2fc) <- samples
print("log2fc germline DEG dataframe!")
head(log2fc)


###3) hierarchical clustering
hclust_matrix <- log2fc %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


rownames(hclust_matrix) <- rownames(log2fc)
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")


library(ComplexHeatmap)
#number of clusters
clustnumb <- 5
library(circlize)
col_fun = colorRamp2(range(hclust_matrix), hcl_palette = "Reds", reverse = TRUE)
p <- Heatmap(hclust_matrix, show_row_names = TRUE, cluster_columns = FALSE, split = clustnumb, heatmap_legend_param = list(title = "Z-score log2FC"), col = col_fun)
#show_row_names = FALSE,

#save the heatmap
pdf(file = snakemake@output[[1]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 8.5) # The height of the plot in inches
draw(p)
dev.off()

#get the identity and number of genes that are in each cluster 
p <- Heatmap(hclust_matrix, show_row_names = FALSE, cluster_columns = FALSE, split = clustnumb)
ht <- draw(p)

clustergenes <- data.frame()

for (u in 1:clustnumb){
	genie <- t(t(row.names(hclust_matrix[row_order(ht)[[u]],])))
	clustdf <- data.frame(SYMBOL = genie, Cluster = rep(u, length(genie)))
	clustergenes <- rbind(clustergenes, clustdf)
}
clustergenes <- merge(clustergenes, XXvsXY, by = "SYMBOL", drop = FALSE )


write.table(clustergenes[order(clustergenes$Cluster),], snakemake@output[[2]], sep = ",", row.names = FALSE )

clustercount <- data.frame()
for (k in 1:clustnumb){
	oneclust <- subset(clustergenes, Cluster == k)
	countclust <- data.frame(Cluster = k, Gene_Number = nrow(oneclust))
	clustercount <- rbind(clustercount, countclust)
}

print(clustercount)



# #make a bar graph of the number of genes in each cluster
library("ggpubr")
pclust <- ggbarplot(clustercount, "Cluster", "Gene_Number",
  color = "Cluster", fill ="Cluster", ylab = "# of shared germline DEGs", xlab = "Cluster #",
  label = TRUE, lab.col = "black", lab.pos = "out", title = "# of genes per cluster") 

pclust <- ggpar(pclust, xticks.by = 1, ylim = c(0,20))
print("making the bar graph")
ggsave(snakemake@output[[3]], plot = pclust + rremove("legend"), width = 4, height = 3)

