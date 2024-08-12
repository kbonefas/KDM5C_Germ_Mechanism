#24.08.11 - To plot the expression of primordial germ cell genes in WT and 5CKO EpiLCs #

###1) Read in the TPM file

#function to format the TPM files
#name is the name in snakemake file
EpiLC_TPM <- read.csv(file = snakemake@input[[1]], sep ="\t", row.names = 1)
genes <-(row.names(EpiLC_TPM))
#remove variant numbers from ensembl gene notation
ENSEMBL<- read.csv(text = genes, sep=".", header = FALSE, col.names = c("ENSEMBL", "var"))
   #change the rownames of the counts file to the ensembl genes without the variant numbers
row.names(EpiLC_TPM) <- ENSEMBL$ENSEMBL

#format the column names by removing extraneous info
colnames(EpiLC_TPM) <- gsub(".TPM", "", colnames(EpiLC_TPM))
head(EpiLC_TPM)




###2) read in the sample info
EpiLC_SI <- read.csv(snakemake@input[[2]], sep = ",")
#get just the male data
#XYEpiLC <- subset(EpiLC_TPM, select =  grep("XY", colnames(EpiLC_TPM), value = TRUE))
#XYEpiLC_SI <- subset(EpiLC_SI, Sex == "XY")
EpiLC_SI$genosex <- paste0(EpiLC_SI$Sex, EpiLC_SI$Genotype)
EpiLC_SI$Sample <- colnames(EpiLC_TPM)



###3) Make a dataframe with the gene info and the ENSEMBL ID
genelist <- data.frame(ENSEMBL = c("ENSMUSG00000024206"), Symbol = c("Rfx2"))

###4) Function to make boxplot of expression in WT and 5cKO EpiLCs
library(ggplot2)
library("ggpubr")
source('code/utilities/colorpalettes.R') #load the colorpalette

#### plot just the male samples, using facet
#generates a tpm plot based on a gene dataframe that has the ensembl IDs in the first column
#ymax - maximum y value
#get the genes that match by ensembl ID
EpiLC_pgc_TPM <- subset(EpiLC_TPM, rownames(EpiLC_TPM) %in% genelist[,1])
#make a new dataframe with  columns being the sample, genotype, and TPM
t_EpiLC_pgc_TPM <- t(EpiLC_pgc_TPM)
print(head(t_EpiLC_pgc_TPM))

plotdf <- data.frame()
for (i in genelist[,1]){
	geneID <- subset(genelist, ENSEMBL == i)
	tempdf <- data.frame(Sample = rownames(t_EpiLC_pgc_TPM), TPM = t_EpiLC_pgc_TPM[,i], ENSEMBL = rep(i, length(rownames(t_EpiLC_pgc_TPM))), Symbol = rep(geneID[,2],  length(rownames(t_EpiLC_pgc_TPM))))

	plotdf <- rbind(plotdf, tempdf)
}


plotdf <- merge(plotdf, EpiLC_SI, by = "Sample")

	#subset for just males
plotdf <- subset(plotdf, plotdf$Sex == "XY")
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
q <- ggpar(q, x.text.angle = 25, font.main = "bold.italic")
	

ggsave(snakemake@output[[1]], plot = q, width = 3, height = 3.5)






