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

#subset for only cells exposed to VA or ESCs (normal culture condition)
ESCEpi_SI_VA <- subset(ESCEpi_SI, ESCEpi_SI$VA != "NO")

#make a variable with the cell type
ESCEpi_SI_VA$Tissue <- ifelse(ESCEpi_SI_VA$Timepoint == "0", "nESC", ifelse(ESCEpi_SI_VA$Timepoint == "48", "EpiLC", ifelse(ESCEpi_SI_VA$Timepoint == "96", "exEpiLC", "uh oh" )))
ESCEpi_SI_VA$GenoTissue <- paste0(ESCEpi_SI_VA$Tissue, "_", ESCEpi_SI_VA$Genotype)

#keep only the needed columns
ESCEpi_SI_VA$Sample <- ESCEpi_SI_VA$ID
ESCEpi_SI_VA <- subset(ESCEpi_SI_VA, select = c("Sample", "GenoTissue"))

print("ESCEpi_SI_VA")
print(ESCEpi_SI_VA)

# ################ AMY and HIP #############
#amygdala
AMY_TPM <- read.csv(file = snakemake@input[["AMY"]], sep ="\t", row.names = 1)
AMY_TPM <- format_ensemb(AMY_TPM)

#remove extraneous info from sample names
colnames(AMY_TPM) <- gsub(".TPM", "", colnames(AMY_TPM))
colnames(AMY_TPM) <- gsub("Sample", "", colnames(AMY_TPM))
print("AMY_TPM")
head(AMY_TPM)



#hippocampus
HIP_TPM <- read.csv(file = snakemake@input[["HIP"]], sep ="\t", row.names = 1)
HIP_TPM <- format_ensemb(HIP_TPM)
#remove extraneous info from sample names
colnames(HIP_TPM) <- gsub(".TPM", "", colnames(HIP_TPM))
colnames(HIP_TPM) <- gsub("Sample", "", colnames(HIP_TPM))
print("HIP_TPM")
head(HIP_TPM)

#make the sample naming scheme in the sample information sheet match the column names of the TPM file
AMYHIP_SI <- read.csv(snakemake@input[["AMYHIP_SI"]], sep =",") 

AMYHIP_SI$Sample <- gsub("_WT", "", AMYHIP_SI$Sample)
AMYHIP_SI$Sample <- gsub("_5cKO", "", AMYHIP_SI$Sample)
print("AMYHIP_SI")
head(AMYHIP_SI)

#add in the amygdala sample information from the ear tag #s
AMY_info <- data.frame(Sample = c("AMY2572", "AMY2622", "AMY2816", "AMY2818", "AMY2879", "AMY2881", "AMY2886", "AMYNT3"), Genotype = c("WT", "5CKO", "5CKO", "5CKO", "WT", "WT", "WT", "5CKO")) 
AMY_info$Reigon <- rep("amygdala", nrow(AMY_info))
AMY_info$Type <- rep("paired.end", nrow(AMY_info))
print("AMY_info")
print(AMY_info)

AMYHIP_SI <- rbind(AMYHIP_SI, AMY_info)
AMYHIP_SI$Tissue <- ifelse(AMYHIP_SI$Reigon == "hippocampus", "HIP", ifelse(AMYHIP_SI$Reigon == "amygdala", "AMY", "uh oh"))
print("AMYHIP_SI new amy ids")
print(AMYHIP_SI)

# rename 5cKO to 5CKO
AMYHIP_SI$Genotype[AMYHIP_SI$Genotype == "5cKO"] <- "5CKO"


AMYHIP_SI$GenoTissue <-  paste0(AMYHIP_SI$Tissue, "_", AMYHIP_SI$Genotype)
AMYHIP_SI <- subset(AMYHIP_SI, select = c("Sample", "GenoTissue"))
print("AMYHIP_SI genotissue")
print(AMYHIP_SI)
# rownames(AMYHIP_SI) <- AMYHIP_SI$ID







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
	plotdf2 <- merge(plotdf, si, by = "Sample")

	return(plotdf2)

}


plotallTPM <- function(genelist){

	#make the plotting df for ESC/EpiLC, AMY, and HIP
	plotdf_ESCEpi <- makeplotdf(genelist, ESCEpi_TPM, ESCEpi_SI_VA)
	print("plotdf_ESCEpi")
	print(head(plotdf_ESCEpi))

	plotdf_AMY <- makeplotdf(genelist, AMY_TPM, AMYHIP_SI)
	print("plotdf_AMY")
	print(head(plotdf_AMY))

	plotdf_HIP <- makeplotdf(genelist, HIP_TPM, AMYHIP_SI)
	print("plotdf_HIP")
	print(head(plotdf_HIP))

	#merge all the dfs together
	plotdf <- rbind(plotdf_ESCEpi, plotdf_AMY, plotdf_HIP)
	
	#order the plotting dataframe
	plotdf$GenoTissue <- factor(plotdf$GenoTissue, levels = c("nESC_WT", "nESC_5CKO", "EpiLC_WT", "EpiLC_5CKO", "exEpiLC_WT", "exEpiLC_5CKO", "AMY_WT", "AMY_5CKO","HIP_WT", "HIP_5CKO"))

	#order genes
	plotdf$Symbol <- factor(plotdf$Symbol, levels = genelist[,2])

	print("All plotting df:")
	print(head(plotdf))

	my_comparisons <- list(c("nESC_WT", "nESC_5CKO"), c("EpiLC_WT", "EpiLC_5CKO"), c("exEpiLC_WT", "exEpiLC_5CKO"), c("AMY_WT", "AMY_5CKO"), c("HIP_WT", "HIP_5CKO"))
	q <- ggboxplot(plotdf, x = 'GenoTissue', y = 'TPM', color = "black", add.params = list(dotsize = 0.1), fill="GenoTissue", 
		add =  "dotplot", xlab = " ", palette = wtKO_pallete) +
    	rremove("legend") +
    	stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.format") 
	q <- ggpar(q, x.text.angle = 50, font.main = "bold.italic")
	
	return(q)
}




#plot all pgc genes of interest:
pgcplot <- plotallTPM(PGCgenes)

ggsave(snakemake@output[[1]], plot = facet(pgcplot, facet.by = "Symbol", nrow = 1), width = 20, height = 3.5)


#plot just a few germline drivers and 2-cell state drivers
PGCgenes_small <- subset(PGCgenes, Symbol %in% c("Dazl", "Stra8",  "Stella (Dppa3)"))
twocellgenes <- data.frame(ENSEMBL = c("ENSMUSG00000075046", "ENSMUSG00000054272", "ENSMUSG00000090714"), Symbol = c("Dux (Duxf3)", "Zscan4c", "Zscan4d"))


pgc_small <- facet(plotallTPM(PGCgenes_small), facet.by = "Symbol", nrow = 1, scales = "free_y")
twocell <- facet(plotallTPM(twocellgenes), facet.by = "Symbol", nrow = 1, scales = "free_y")


library("gridExtra")
ggsave(snakemake@output[[2]], plot = grid.arrange(grobs = list(pgc_small, twocell), nrow = 2), width = 12, height = 8)



# #piRNA genes
# piRNAgenes <- data.frame(ENSEMBL = c("ENSMUSG00000021758", "ENSMUSG00000033644", "ENSMUSG00000009628", "ENSMUSG00000035517"), Symbol = c("Ddx4 (Mvh)", "Piwil2 (Mili)", "Tex15", "Tdrd7"))

# piRNA <- facet(plotEpiLCTPM(piRNAgenes, 10), facet.by = "Symbol", nrow = 1)
# ggsave(snakemake@output[[3]], plot = piRNA, width = 5, height = 2.5)







