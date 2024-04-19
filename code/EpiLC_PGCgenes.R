# To plot the expression of primordial germ cell genes in WT and 5CKO EpiLCs #

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
PGCgenes <- data.frame(ENSEMBL = c("ENSMUSG00000010592","ENSMUSG00000029848", "ENSMUSG00000021758", "ENSMUSG00000046323","ENSMUSG00000025492"), Symbol = c("Dazl", "Stra8",  "Mvh (Ddx4)", "Stella (Dppa3)", "Fragilis (Ifitm3)"))
#ENSMUSG00000038151 - "Prdm1 (Blimp1)" (not a DEG)
#ENSMUSG00000005672 - "c-Kit" (only XX DEG)

###4) Function to make boxplot of expression in WT and 5cKO EpiLCs
library(ggplot2)
library("ggpubr")
source('code/utilities/colorpalettes.R') #load the colorpalette

# tpmboxplot <- function(df, gene){
# 		my_comparisons <- list(c("XYWT", "XY5cKO"), c("XXWT", "XX5cKO"), c("XXWT", "XX5cHET"))

# 		pl <- ggboxplot(df, x = 'genosex', y = 'TPM', color = "black", fill="genosex", 
#     	title = gene, 
#     	add =  "dotplot", add.params = list(size = 0.7), xlab = " ", 
#     	palette = EpiLCpalette) +
#     	rremove("legend") + scale_x_discrete(labels = c("XY 5cKO", "XY WT", "XX WT", "XX 5cHET", "XX 5cKO")) +
#         stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.format") 
# 		pl <- ggpar(pl, x.text.angle = 25, font.main = "bold.italic")
	
# }



# ###4) Make dataframe and plot with the gene of interest
# 	#for every gene, make a df

# #an empty list to store the plots in
# TPMplot <- list()
# count <- 1

# for (i in unique(PGCgenes$Symbol)){
# 	ens <- subset(PGCgenes, Symbol == i)
# 	#get the genes that match by ensembl ID
# 	GOI <- subset(EpiLC_TPM, rownames(EpiLC_TPM) == ens[,1])
# 	print(head(GOI))

# 	#make a new dataframe with  columns being the sample, genotype, and TPM
# 	tGOI <- t(GOI)
# 	plotdf <- data.frame(Sample = rownames(tGOI), TPM = tGOI[,1])
# 	plotdf <- merge(plotdf, EpiLC_SI, by = "Sample")
# 	#order the factor levels so WT plots first
# 	plotdf$genosex <- factor(plotdf$genosex, levels = names(EpiLCpalette))
# 	print(paste0("Dataframe for ", i))
# 	print(head(plotdf))


# 	#generate a plot for the gene
# 	q <- tpmboxplot(plotdf, i)
# 	TPMplot[[count]] <- q
# 	count <- count + 1

# }

# library("gridExtra")
# ggsave(snakemake@output[[1]], plot = grid.arrange(grobs = TPMplot, nrow = 2), width = 9, height = 9)

#### plot just the male samples, using facet

#get the genes that match by ensembl ID
EpiLC_pgc_TPM <- subset(EpiLC_TPM, rownames(EpiLC_TPM) %in% PGCgenes[,1])
#make a new dataframe with  columns being the sample, genotype, and TPM
t_EpiLC_pgc_TPM <- t(EpiLC_pgc_TPM)
print(head(t_EpiLC_pgc_TPM))

plotdf <- data.frame()

for (i in PGCgenes[,1]){
	geneID <- subset(PGCgenes, ENSEMBL == i)
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
plotdf$Symbol <- factor(plotdf$Symbol, levels = PGCgenes[,2])

print("All germ genes df")
print(head(plotdf))

my_comparisons <- list(c("WT", "5CKO"))
q <- ggboxplot(plotdf, x = 'Genotype', y = 'TPM', color = "black", fill="Genotype", 
	add.params = list(size = 1.25), add =  "dotplot", xlab = " ", palette = EpiLC_XY_palette) +
    rremove("legend") +
    stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.format") 
q <- ggpar(q, x.text.angle = 25, ylim = c(0,25), font.main = "bold.italic")
	

ggsave(snakemake@output[[1]], plot = facet(q, facet.by = "Symbol", nrow = 1), width = 6, height = 2.5)










