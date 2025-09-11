# To plot the expression of primordial germ cell genes in WT and 5CKO ESCs, EpiLCs, and brain #

#need a dataframe with the columns gene column, tissue, genotype, TPM, tissgeno
	#then subset based on the gene and make a boxplot

#read in the TPM file
#remove variant numbers on gene names
#get which sample is which - read in the sample info, merge transposed df with the sample names


#remove transcript variants from ensembl gene names (rownames of df)
format_ensemb <- function(df){
	genes <-(row.names(df))
	#remove variant numbers from ensembl gene notation
	ENSEMBL<- read.csv(text = genes, sep=".", header = FALSE, col.names = c("ENSEMBL", "var"))
	#change the rownames of the counts file to the ensembl genes without the variant numbers
	df2 <- df
	row.names(df2) <- ENSEMBL$ENSEMBL
	return(df2)

}

################ ESC to exEpiLC #############
#1) Read in the TPM file

#function to format the TPM files
#name is the name in snakemake file
ESCEpi_TPM <- read.csv(file = snakemake@input[["ESCEpi"]], sep ="\t", row.names = 1)
ESCEpi_TPM <- ESCEpi_TPM[4:ncol(ESCEpi_TPM)] 
print("ESCEpi_TPM")
print(head(ESCEpi_TPM))

####2)get the sample information

#format the TPM data
#reorder the TPM based on the samples
ESCEpi_SI <- read.csv(snakemake@input[["ESCEpi_SI"]], sep =",") 
rownames(ESCEpi_SI) <- ESCEpi_SI$ID

#subset for only cells exposed to VA (normal culture condition)
ESCEpi_SI_VA <- subset(ESCEpi_SI, ESCEpi_SI$VA == "VA")

#make a variable with the cell type
ESCEpi_SI_VA$Tissue <- ifelse(ESCEpi_SI_VA$Timepoint == "0", "nESC", ifelse(ESCEpi_SI_VA$Timepoint == "48", "EpiLC", ifelse(ESCEpi_SI_VA$Timepoint == "96", "exEpiLC", "uh oh" )))
ESCEpi_SI_VA$GenoTissue <- paste0(ESCEpi_SI_VA$Tissue, "_", ESCEpi_SI_VA$Genotype)

#keep only the needed columns
ESCEpi_SI_VA <- subset(ESCEpi_SI_VA, select = c("Sample", "GenoTissue"))

print("ESCEpi_SI_VA")
print(ESCEpi_SI_VA)

################ AMY and HIP #############
AMY_TPM <- read.csv(file = snakemake@input[["AMY"]], sep ="\t", row.names = 1)
ESCEpi_TPM <- ESCEpi_TPM[4:ncol(ESCEpi_TPM)] 
print("ESCEpi_TPM")
print(head(ESCEpi_TPM))





###3) Make a dataframe with the gene info and the ENSEMBL ID
PGCgenes <- data.frame(ENSEMBL = c("ENSMUSG00000010592","ENSMUSG00000029848", "ENSMUSG00000021758", "ENSMUSG00000046323","ENSMUSG00000025492"), Symbol = c("Dazl", "Stra8",  "Mvh (Ddx4)", "Stella (Dppa3)", "Fragilis (Ifitm3)"))


###4) Function to make boxplot of expression in WT and 5cKO
library(ggplot2)
library("ggpubr")
source('code/utilities/colorpalettes.R') #load the colorpalette


#list of genes, tpm dataframe, and sample info
makeplotdf <- function(genelist, tpm, si){
	#get the genes that match by ensembl ID
	tpm_genes <- subset(tpm, rownames(tpm) %in% genelist[,1])
	#make a new dataframe with  columns being the sample, genotype, and TPM, and tissue
	t_tpm_genes <- t(tpm_genes)
	print("t_tpm_genes")
	print(head(t_tpm_genes))

	plotdf <- data.frame()

	for (i in genelist[,1]){
		geneID <- subset(genelist, ENSEMBL == i)
		tempdf <- data.frame(Sample = rownames(t_tpm_genes), TPM = t_tpm_genes[,i], ENSEMBL = rep(i, length(rownames(t_tpm_genes))), Symbol = rep(geneID[,2],  length(rownames(t_tpm_genes))))

		plotdf <- rbind(plotdf, tempdf)
	}
	plotdf <- merge(plotdf, si, by = "Sample")
	return(plotdf)

}





plotallTPM <- function(genelist, ymax){


	#make the plotting df for ESC/EpiLC, AMY, and HIP
	plotdf_ESCEpi <- makeplotdf(genelist, ESCEpi_TPM, ESCEpi_SI_VA)
	print(head(plotdf_ESCEpi))





	#order the factor levels so WT plots first and rename 5cKO
	plotdf$Genotype[plotdf$Genotype == "5cKO"] <- "5CKO"
	plotdf$Genotype <- factor(plotdf$Genotype, levels = c("WT", "5CKO"))

	#order genes
	plotdf$Symbol <- factor(plotdf$Symbol, levels = genelist[,2])

	print("All plotting df:")
	print(head(plotdf))

	my_comparisons <- list(c("WT", "5CKO"))
	q <- ggboxplot(plotdf, x = 'Genotype', y = 'TPM', color = "black", add.params = list(size = 1.25), fill="Genotype", 
		add =  "dotplot", xlab = " ", palette = EpiLC_XY_palette) +
    	rremove("legend") +
    	stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.format") 
	q <- ggpar(q, x.text.angle = 25, ylim = c(0,ymax), font.main = "bold.italic")
	
	return(q)
}






#### plot just the male samples, using facet
#generates a tpm plot based on a gene dataframe that has the ensembl IDs in the first column
#ymax - maximum y value
# plotEpiLCTPM <- function(genelist, ymax){
# 	#get the genes that match by ensembl ID
# 	EpiLC_pgc_TPM <- subset(EpiLC_TPM, rownames(EpiLC_TPM) %in% genelist[,1])
# 	#make a new dataframe with  columns being the sample, genotype, and TPM
# 	t_EpiLC_pgc_TPM <- t(EpiLC_pgc_TPM)
# 	print(head(t_EpiLC_pgc_TPM))

# 	plotdf <- data.frame()

# 	for (i in genelist[,1]){
# 		geneID <- subset(genelist, ENSEMBL == i)
# 		tempdf <- data.frame(Sample = rownames(t_EpiLC_pgc_TPM), TPM = t_EpiLC_pgc_TPM[,i], ENSEMBL = rep(i, length(rownames(t_EpiLC_pgc_TPM))), Symbol = rep(geneID[,2],  length(rownames(t_EpiLC_pgc_TPM))))

# 		plotdf <- rbind(plotdf, tempdf)
# 	}


# 	plotdf <- merge(plotdf, EpiLC_SI, by = "Sample")

# 	#subset for just males
# 	plotdf <- subset(plotdf, plotdf$Sex == "XY")
# 	#order the factor levels so WT plots first and rename 5cKO
# 	plotdf$Genotype[plotdf$Genotype == "5cKO"] <- "5CKO"
# 	plotdf$Genotype <- factor(plotdf$Genotype, levels = c("WT", "5CKO"))

# 	#order genes
# 	plotdf$Symbol <- factor(plotdf$Symbol, levels = genelist[,2])

# 	print("All plotting df:")
# 	print(head(plotdf))

# 	my_comparisons <- list(c("WT", "5CKO"))
# 	q <- ggboxplot(plotdf, x = 'Genotype', y = 'TPM', color = "black", add.params = list(size = 1.25), fill="Genotype", 
# 		add =  "dotplot", xlab = " ", palette = EpiLC_XY_palette) +
#     	rremove("legend") +
#     	stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.format") 
# 	q <- ggpar(q, x.text.angle = 25, ylim = c(0,ymax), font.main = "bold.italic")
	
# 	return(q)
# }



#plot all pgc genes of interest:
pgcplot <- plotEpiLCTPM(PGCgenes, 25)

ggsave(snakemake@output[[1]], plot = facet(pgcplot, facet.by = "Symbol", nrow = 1), width = 6, height = 2.5)


#plot just a few germline drivers and 2-cell state drivers
PGCgenes_small <- subset(PGCgenes, Symbol %in% c("Dazl", "Stra8",  "Stella (Dppa3)"))
twocellgenes <- data.frame(ENSEMBL = c("ENSMUSG00000075046", "ENSMUSG00000054272", "ENSMUSG00000090714"), Symbol = c("Dux (Duxf3)", "Zscan4c", "Zscan4d"))


pgc_small <- facet(plotEpiLCTPM(PGCgenes_small, 25), facet.by = "Symbol", nrow = 1)
twocell <- facet(plotEpiLCTPM(twocellgenes, 10), facet.by = "Symbol", nrow = 1)


library("gridExtra")
ggsave(snakemake@output[[2]], plot = grid.arrange(grobs = list(pgc_small, twocell), nrow = 1), width = 7, height = 2.5)



#piRNA genes
piRNAgenes <- data.frame(ENSEMBL = c("ENSMUSG00000021758", "ENSMUSG00000033644", "ENSMUSG00000009628", "ENSMUSG00000035517"), Symbol = c("Ddx4 (Mvh)", "Piwil2 (Mili)", "Tex15", "Tdrd7"))

piRNA <- facet(plotEpiLCTPM(piRNAgenes, 10), facet.by = "Symbol", nrow = 1)
ggsave(snakemake@output[[3]], plot = piRNA, width = 5, height = 2.5)







