#2025.05.22 - heatmap of all germline genes in WT and 5CKO cells as they differentiate from ESC to exEpiLCs

#read in germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")
print(head(germ))

#read in TPM
tpm_all <- read.csv(snakemake@input[[2]], sep ="\t", row.names = 1)
tpm <- tpm_all[4:ncol(tpm_all)] 
print(head(tpm))

#format the TPM data
#reorder the TPM based on the samples
SampleInfo <- read.csv(snakemake@input[[3]], sep =",") 
print(SampleInfo)
rownames(SampleInfo) <- SampleInfo$ID

#make a group variable that has the genotype, time point, and RA status
SampleInfo$group <- paste0(SampleInfo$Genotype,"_", SampleInfo$Timepoint, SampleInfo$VA)

#get average TPM
avgtpm_df <- list()
count = 1
for (g in unique(SampleInfo$group)){
	#get columns with samples of that genotype/RA/Time point
	samp <- subset(SampleInfo, group == g)
	samp <- samp$ID

	onegeno <- subset(tpm, select = samp)
	# print("one geno")
	# print(head(onegeno))
	#get the average TPM
	avgTPM <- data.frame(row.names = rownames(tpm), TPM = rowMeans(onegeno))
	names(avgTPM)[names(avgTPM) == 'TPM'] <- g 
	avgtpm_df[[count]] <- avgTPM
	count = count + 1
}

library(dplyr)
avgtpm_df <- bind_cols(avgtpm_df)
print("avg TPM dataframe!")
# head(avgtpm_df)

#subset TPM for germline genes
germ_avg <- subset(avgtpm_df, rownames(avgtpm_df) %in% germ$ENSEMBL)

#subset for WT ESC values < 1 (unexpressed in ESCs)
germ_avg <- subset(germ_avg, germ_avg$WT_0 < 1)


#remove genes unexpressed in any sample (> or equal to 1 in at least one column)
germ_avg <- germ_avg[rowSums(germ_avg >= 1) >= 1, ]

#order the columns
germ_avg <- germ_avg[, c("WT_0", "5CKO_0", "WT_48NO", "5CKO_48NO", "WT_48VA", "5CKO_48VA", "WT_96NO", "5CKO_96NO", "WT_96VA", "5CKO_96VA")] 
print(head(germ_avg))


#transform TPM into log2(TPM+1)
germ_log2 <- log2(germ_avg+1)
print("log2+1")
print(head(germ_log2))



#plot heatmap of germline genes

library(ComplexHeatmap)
#number of clusters

#maximum value 
print(paste("max value", max(germ_log2)))
library(circlize)
col_fun = colorRamp2(0:max(germ_log2), hcl_palette = "Reds", reverse = TRUE)


p <- Heatmap(germ_log2, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, heatmap_legend_param = list(title = "log2 TPM + 1", legend_direction = "horizontal"), col = col_fun, column_split  = c(rep("0", 2), rep("48 No Vit A", 2), rep("48 Vit A", 2), rep("96 No Vit A", 2), rep("96 Vit A", 2)), column_title = "germline gene expression")


pdf(file = snakemake@output[[1]], width = 8, height = 10)
	draw(p)
dev.off()



#another heat map that groups WT time points together and 5CKO together
#change column order
germ_log2 <- germ_log2[, c("WT_0", "WT_48NO", "WT_48VA", "WT_96NO", "WT_96VA", "5CKO_0", "5CKO_48NO", "5CKO_48VA", "5CKO_96NO", "5CKO_96VA")] 

p <- Heatmap(germ_log2, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, heatmap_legend_param = list(title = "log2 TPM + 1", legend_direction = "horizontal"), col = col_fun, column_split  = c(rep("WT", 5), rep("KO", 5)), column_title = "germline gene expression")


pdf(file = snakemake@output[[2]], width = 8, height = 10)
	draw(p)
dev.off()



#### plot z-scores instead

#germ_avg - average tpm of all germline genes

# get scaled FPKM 
library(tidyverse)
print("germ_avg")
print(head(germ_avg))

#scaled matrix
hclust_matrix <- as.matrix(germ_avg) %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

# rownames(hclust_matrix) <- germ_avg$gene_short_name
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
print("hclust_matrix")
print(head(hclust_matrix))

#plot heatmap of scaled TPM

clustnumb <- 4
library(circlize)
col_fun = colorRamp2(0:max(hclust_matrix), hcl_palette = "Reds", reverse = TRUE)

# #print(head(hclust_matrix_t))
# p <- Heatmap(hclust_matrix, column_split = clustnumb, heatmap_legend_param = list(title = "Z-score FPKM"), col = col_fun)


p <- Heatmap(hclust_matrix, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, heatmap_legend_param = list(title = "Z-score TPM", legend_direction = "horizontal"), col = col_fun, column_split  = c(rep("0", 2), rep("48 No Vit A", 2), rep("48 Vit A", 2), rep("96 No Vit A", 2), rep("96 Vit A", 2)), column_title = "germline gene expression")


pdf(file = snakemake@output[[3]], width = 8, height = 10)
	draw(p)
dev.off()


###zscore
#another heat map that groups WT time points together and 5CKO together
#change column order
hclust_matrix <- hclust_matrix[, c("WT_0", "WT_48NO", "WT_48VA", "WT_96NO", "WT_96VA", "5CKO_0", "5CKO_48NO", "5CKO_48VA", "5CKO_96NO", "5CKO_96VA")] 

p <- Heatmap(hclust_matrix, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, heatmap_legend_param = list(title = "Z-score TPM", legend_direction = "horizontal"), col = col_fun, column_split  = c(rep("WT", 5), rep("KO", 5)), column_title = "germline gene expression")


#save the heatmap
pdf(file = snakemake@output[[4]],   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6) # The height of the plot in inches
draw(p)
dev.off() 

