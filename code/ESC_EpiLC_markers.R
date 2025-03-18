#see how the expression of germline DEGs is affected by RA
tpm_all <- read.csv(snakemake@input[[1]], sep ="\t", row.names = 1)
tpm <- tpm_all[4:ncol(tpm_all)] 
print(head(tpm))

#format the TPM data
#reorder the TPM based on the samples
SampleInfo <- read.csv(snakemake@input[[2]], sep =",") 
rownames(SampleInfo) <- SampleInfo$ID

#make a group variable that has the genotype, time point, and RA status
SampleInfo$group <- paste0(SampleInfo$Genotype,"_", SampleInfo$Timepoint, SampleInfo$RA)

#for each group, calculate the average tpm
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
head(avgtpm_df)


#list of marker genes to look for
#should be high in ESCs, low in WT EpiLCs
naivegenes <- data.frame(Symbol = c("Esrrb", "Fgf4", "Klf2", "Klf4", "Klf5", "Nr0b1", "Nr5a2", "Prdm14", "Sox13", "Tbx3", "Tcl1", "Tet2", "Tfcp2l1", "Zfp42"), ENSEMBL = c("ENSMUSG00000021255", "ENSMUSG00000050917", "ENSMUSG00000055148", "ENSMUSG00000003032", "ENSMUSG00000005148", "ENSMUSG00000025056", "ENSMUSG00000026398", "ENSMUSG00000042414", "ENSMUSG00000070643", "ENSMUSG00000018604", "ENSMUSG00000041359", "ENSMUSG00000040943", "ENSMUSG00000026380", "ENSMUSG00000051176"))

#should be low in ESCs, high in WT EpiLC
primedgenes <- data.frame(Symbol = c("Cdh2", "Dnmt3a", "Dnmt3b", "Etv4", "Etv5", "Foxd3", "Lef1", "Otx2", "Pou3f1", "Sox3", "Tcf15", "Fgf5"), ENSEMBL = c("ENSMUSG00000024304", "ENSMUSG00000020661", "ENSMUSG00000027478", "ENSMUSG00000017724", "ENSMUSG00000013089", "ENSMUSG00000067261", "ENSMUSG00000027985", "ENSMUSG00000021848", "ENSMUSG00000090125", "ENSMUSG00000045179", "ENSMUSG00000068079", "ENSMUSG00000029337"))

genelist <- list(naivegenes, primedgenes)
names(genelist) <- c("naive", "primed")


heat <- list()

for (i in 1:length(genelist)){

	df <- genelist[[i]]
	#get the tpm for the genes in each list
	tpm_GOI <- subset(avgtpm_df, rownames(avgtpm_df) %in% df$ENSEMBL)
	#replace the ensembl names with the symbol names
	row.names(df) <- df$ENSEMBL
	tpm_GOI <- merge(tpm_GOI, df, by = "row.names")

	row.names(tpm_GOI) <- tpm_GOI$Symbol
	tpm_GOI <- subset(tpm_GOI, select = c(-Symbol, -ENSEMBL))

	#order the columns
	tpm_GOI <- tpm_GOI[, c("WT_0", "5CKO_0", "WT_48NO", "5CKO_48NO", "WT_48RA", "5CKO_48RA", "WT_96NO", "5CKO_96NO", "WT_96RA", "5CKO_96RA")] 
	print(head(tpm_GOI))

	tpm_GOI <- log2(tpm_GOI+1)

	# hclust_matrix <- tpm_GOI %>% 
 	#  	# transpose the matrix so genes are as columns
  	# t() %>% 
 	# 	 # apply scaling to each column of the matrix (genes)
  	# scale() %>% 
  	# 	# transpose back so genes are as 	rows again
  	# t()



	library(ComplexHeatmap)
	#number of clusters
	library(circlize)
	col_fun = colorRamp2(0:10, hcl_palette = "Reds", reverse = TRUE)


	p <- Heatmap(tpm_GOI, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title = "log2 TPM + 1", legend_direction = "horizontal"), col = col_fun, column_split  = c(rep("0", 2), rep("48 No Vit A", 2), rep("48 Vit A", 2), rep("96 No Vit A", 2), rep("96 Vit A", 2)), column_title = paste(names(genelist)[i], "pluripotency genes"))

	heat[[i]] <- draw(p)

	pdf(file = snakemake@output[[i]], width = 8, height = 4) # 
		draw(p)
	dev.off()

}

#1) change gene names
#2) Fix scaling
#3) Fix order

# library("gridExtra")
# library("ggpubr")
# ggsave(snakemake@output[[1]], plot = grid.arrange(grobs = heat, nrow = length(genelist)), width = 10, height = 8.5)







