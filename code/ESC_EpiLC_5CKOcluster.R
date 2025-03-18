#23.10.10 Visualizing germline genes dysregulated with Kdm5c loss in ESCs and EpiLCs

##1) Read in the ESC and EpiLC germline DEGs
germ <- read.csv(snakemake@input[[1]], sep = ",") 

samples <- c("ESC_0", "EpiLC_48NO", "EpiLC_48RA", "EpiLC_96NO", "EpiLC_96RA")
restabs <- list()
germDEGs <- list()
allgermDEGs <- data.frame()
for (s in 1:length(samples)){
	#read in the dataframe
	df <- read.csv(snakemake@input[[s+1]], sep = ",", row.names = 1)
	df$ENSEMBL <- rownames(df)
	#get the germ genes
	germDEG <- subset(df, df$ENSEMBL %in% germ$ENSEMBL & df$padj < snakemake@params[["alpha"]] & df$log2FoldChange > 0)
	germDEG <- merge(germ, germDEG, by = "ENSEMBL")

	restabs[[s]] <- df
	germDEGs[[s]] <- germDEG
	#add to DEGs for the germline gene list
	allgermDEGs <- rbind(allgermDEGs, germDEG[,1:2])
}


names(germDEGs) <- samples
names(restabs) <- samples

print(head(restabs[["ESC_0"]]))

#remove duplicate genes
allgermDEGs <- unique(allgermDEGs)
print(head(allgermDEGs))

##2) make a dataframe with the # of germline genes for each sample
plotdf <- data.frame(Samples = samples, germline_DEGs = sapply(germDEGs, nrow), Time = c("0", "48", "48", "96", "96"), RA = c("-", "-", "+", "-", "+"))
#order the samples for plotting 
plotdf$Samples <- factor(plotdf$Samples, levels = samples)
print(plotdf)

##4) Plot the total # of germline DEGs in a histogram
source("code/utilities/colorpalettes.R")
library("ggpubr")
p <- ggbarplot(plotdf, "Samples", "germline_DEGs",
  color = "Time", fill ="Time", palette = RAtimes, ylab = "# of germline DEGs", xlab = " ",
  label = TRUE, lab.col = "black", lab.pos = "out", title = "# of Germline DEGs") 
  #scale_fill_manual(values = c("+" = "black", "-" = "white"))


p <- ggpar(p, x.text.angle = 25, ylim = c(0,225))
ggsave(snakemake@output[[1]], plot = p + rremove("legend"), width = 4, height = 4)


#see how the expression of germline DEGs is affected by RA
tpm_all <- read.csv(snakemake@input[[7]], sep ="\t", row.names = 1)
tpm <- tpm_all[4:ncol(tpm_all)] 
print(head(tpm))
#subset tpm for just the germline DEGs
tpm_germ <- subset(tpm, rownames(tpm) %in% allgermDEGs$ENSEMBL)

#reorder the TPM based on the samples
SampleInfo <- read.csv(snakemake@input[[8]], sep =",") 
rownames(SampleInfo) <- SampleInfo$ID

#make a group variable that has the genotype, time point, and RA status
SampleInfo$group <- paste0(SampleInfo$Genotype,"_", SampleInfo$Timepoint, SampleInfo$RA)

#do just the 5cKOs 
SampleInfo_5C <- subset(SampleInfo, Genotype == "5CKO" )

#for each group, calculate the average tpm

avg5Cgermdegs <- list()
count = 1
for (g in unique(SampleInfo_5C$group)){
	#get columns with samples of that genotype/RA/Time point
	samp <- subset(SampleInfo_5C, group == g)
	samp <- samp$ID

	onegeno <- subset(tpm_germ, select = samp)
	#get the average TPM
	avgTPM <- data.frame(row.names = rownames(tpm_germ), TPM = rowMeans(onegeno))
	names(avgTPM)[names(avgTPM) == 'TPM'] <- g 
	avg5Cgermdegs[[count]] <- avgTPM
	count = count + 1
}

library(dplyr)
avg5Cgermdegs <- bind_cols(avg5Cgermdegs)
print("avg germline DEG dataframe!")
head(avg5Cgermdegs)

hclust_matrix <- avg5Cgermdegs %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

############################################################################################
#make one version transposed for plotting, the other normal for downstream analyses
hclust_matrix_t <- avg5Cgermdegs %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale()

colnames(hclust_matrix_t) <- rownames(avg5Cgermdegs)
rownames(hclust_matrix) <- rownames(avg5Cgermdegs)
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")


library(ComplexHeatmap)
#number of clusters
clustnumb <- 8
library(circlize)
col_fun = colorRamp2(range(hclust_matrix), hcl_palette = "Reds", reverse = TRUE)
#print(head(hclust_matrix_t))
p <- Heatmap(hclust_matrix_t, show_column_names = FALSE, cluster_rows = FALSE, column_split = clustnumb, heatmap_legend_param = list(title = "Z-score TPM", legend_direction = "horizontal"), col = col_fun, row_order = c("5CKO_0", "5CKO_48NO", "5CKO_48RA", "5CKO_96NO", "5CKO_96RA"))

#save the heatmap
pdf(file = snakemake@output[[2]],   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 3) # The height of the plot in inches
draw(p)
dev.off()

#get the identity and number of genes that are in each cluster 
p <- Heatmap(hclust_matrix, show_row_names = FALSE, cluster_columns = FALSE, split = clustnumb)
ht <- draw(p)

clustergenes <- data.frame()

for (u in 1:clustnumb){
	genie <- t(t(row.names(hclust_matrix[row_order(ht)[[u]],])))
	clustdf <- data.frame(ENSEMBL = genie, Cluster = rep(u, length(genie)))
	clustergenes <- rbind(clustergenes, clustdf)
}
clustergenes <- merge(clustergenes, allgermDEGs, by = "ENSEMBL", drop = FALSE )


write.table(clustergenes[order(clustergenes$Cluster),], snakemake@output[[3]], sep = ",", row.names = FALSE )

clustercount <- data.frame()
for (k in 1:clustnumb){
	oneclust <- subset(clustergenes, Cluster == k)
	countclust <- data.frame(Cluster = k, Gene_Number = nrow(oneclust))
	clustercount <- rbind(clustercount, countclust)
}

print(clustercount)



#make a bar graph of the number of genes in each cluster
library("ggpubr")
pclust <- ggbarplot(clustercount, "Cluster", "Gene_Number",
  color = "Cluster", fill ="Cluster", ylab = "# of germline DEGs", xlab = "Cluster #",
  label = TRUE, lab.col = "black", lab.pos = "out", title = "# of genes per cluster") 

pclust <- ggpar(pclust, xticks.by = 1, ylim = c(0,150))
print("making the bar graph TPM")
ggsave(snakemake@output[[4]], plot = pclust + rremove("legend"), width = 3, height = 3)


#get which cluster brain DEGs fall into
germHIP <- read.csv(snakemake@input[[9]], sep = ",")
germAMY <- read.csv(snakemake@input[[10]], sep = ",")
germDEGlist <- list(germHIP, germAMY)


brainclust <- function(df, i){
	brainDEGclust <- merge(df, clustergenes, by = "ENSEMBL")
	write.table(brainDEGclust, snakemake@output[[i+5]], sep = ",", row.names = FALSE )
} 

for (g in 1:length(germDEGlist)){
	brainclust(germDEGlist[[g]], g)
}




###########################
# genes only dysregulated with RA in 96 hrs vs no RA

RA96only <- subset(germDEGs[["EpiLC_96RA"]], !(germDEGs[["EpiLC_96RA"]][, "ENSEMBL"] %in% germDEGs[["EpiLC_96NO"]][, "ENSEMBL"]))
head(RA96only)
nrow(RA96only)

write.table(RA96only, snakemake@output[[5]], sep = ",", row.names = FALSE )


## 2024.07.01 Get the transcript start and end site of each gene in the clusters for bedtools plotting of KDM5C binding
#get list of genes for each cluster then get the TSS and TES
source('code/utilities/GeneTSSandTES.R')
for (i in 1:clustnumb){
	genelist <- subset(clustergenes, Cluster == i)$ENSEMBL
	geneTSSandTES(genelist, i+7)
}


