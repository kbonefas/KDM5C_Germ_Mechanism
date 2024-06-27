### 24.04.19 Assess expression of ESC and EpiLC genes in WT and 5CKO to ensure differentiation occurred properly
	# use genes in this paper as a starting point: https://pubmed.ncbi.nlm.nih.gov/38177678/
	#make heatmap of TPM in WT and 5CKO

###1) list of naive and primed pluripotency genes

#should be low in WT EpiLCs
naivegenes <- data.frame(Symbol = c("Esrrb", "Fgf4", "Klf2", "Klf4", "Klf5", "Nr0b1", "Nr5a2", "Prdm14", "Sox13", "Tbx3", "Tcl1", "Tet2", "Tfcp2l1", "Zfp42"), ENSEMBL = c("ENSMUSG00000021255", "ENSMUSG00000050917", "ENSMUSG00000055148", "ENSMUSG00000003032", "ENSMUSG00000005148", "ENSMUSG00000025056", "ENSMUSG00000026398", "ENSMUSG00000042414", "ENSMUSG00000070643", "ENSMUSG00000018604", "ENSMUSG00000041359", "ENSMUSG00000040943", "ENSMUSG00000026380", "ENSMUSG00000051176"))

#should be high in WT EpiLC
primedgenes <- data.frame(Symbol = c("Cdh2", "Dnmt3a", "Dnmt3b", "Etv4", "Etv5", "Foxd3", "Lef1", "Otx2", "Pou3f1", "Sox3", "Tcf15", "Fgf5"), ENSEMBL = c("ENSMUSG00000024304", "ENSMUSG00000020661", "ENSMUSG00000027478", "ENSMUSG00000017724", "ENSMUSG00000013089", "ENSMUSG00000067261", "ENSMUSG00000027985", "ENSMUSG00000021848", "ENSMUSG00000090125", "ENSMUSG00000045179", "ENSMUSG00000068079", "ENSMUSG00000029337"))

#genes for core pluripotency, not expected to change with ESC to EpiLCs
coregenes <- data.frame(Symbol = c("Lin28a", "Nanog", "Pou5f1", "Sall4", "Sox2", "Utf1", "Zfp281"), ENSEMBL = c("ENSMUSG00000050966", "ENSMUSG00000012396", "ENSMUSG00000024406", "ENSMUSG00000027547", "ENSMUSG00000074637", "ENSMUSG00000047751", "ENSMUSG00000041483"))


###2) Read in the TPM file

#format the TPM files
#name is the name in snakemake file
EpiLC_TPM <- read.csv(file = snakemake@input[[1]], sep ="\t", row.names = 1)
genes <-(row.names(EpiLC_TPM))
#remove variant numbers from ensembl gene notation
ENSEMBL<- read.csv(text = genes, sep=".", header = FALSE, col.names = c("ENSEMBL", "var"))
   #change the rownames of the counts file to the ensembl genes without the variant numbers
row.names(EpiLC_TPM) <- ENSEMBL$ENSEMBL

#format the column names by removing extraneous info
colnames(EpiLC_TPM) <- gsub(".TPM", "", colnames(EpiLC_TPM))


###3) read in the sample info
EpiLC_SI <- read.csv(snakemake@input[[2]], sep = ",")
#get just the male data
EpiLC_TPM <- subset(EpiLC_TPM, select =  grep("XY", colnames(EpiLC_TPM), value = TRUE))
EpiLC_SI <- subset(EpiLC_SI, Sex == "XY")
EpiLC_SI$Sample <- colnames(EpiLC_TPM)
print(head(EpiLC_TPM))



###4) plot the tpm of naive and primed pluripotency genes in WT vs 5CKO
#put all the genes together
naivegenes$type <- rep("naive", nrow(naivegenes))
primedgenes$type <- rep("primed", nrow(primedgenes))

markergenes <- rbind(naivegenes, primedgenes)
print(head(markergenes))


#get the genes that match by ensembl ID
EpiLC_marker_TPM <- subset(EpiLC_TPM, rownames(EpiLC_TPM) %in% markergenes[,2])
#make a new dataframe with  columns being the sample, genotype, and TPM
t_EpiLC_marker_TPM <- t(EpiLC_marker_TPM)
print(head(t_EpiLC_marker_TPM))

plotdf <- data.frame()

#ENSEMBL is the second column
for (i in markergenes[,2]){
	geneID <- subset(markergenes, ENSEMBL == i)
	tempdf <- data.frame(Sample = rownames(t_EpiLC_marker_TPM), TPM = t_EpiLC_marker_TPM[,i], ENSEMBL = rep(i, length(rownames(t_EpiLC_marker_TPM))), Symbol = rep(geneID[,1],  length(rownames(t_EpiLC_marker_TPM))), type = rep(geneID[,3],  length(rownames(t_EpiLC_marker_TPM))))

	plotdf <- rbind(plotdf, tempdf)
}

plotdf <- merge(plotdf, EpiLC_SI, by = "Sample")
print(head(plotdf))
#order the factor levels so WT plots first and rename 5cKO
plotdf$Genotype[plotdf$Genotype == "5cKO"] <- "5CKO"
plotdf$Genotype <- factor(plotdf$Genotype, levels = c("WT", "5CKO"))


library("ggpubr")
source('code/utilities/colorpalettes.R') #load the colorpalette

#make separate plots for the gene type and combine together for final plot
TPMplot <- list()
genetype <- unique(plotdf$type)

for (i in 1:length(genetype)){
	plotdf_2 <- subset(plotdf, type == genetype[i])
	my_comparisons <- list(c("WT", "5CKO"))
	q <- ggboxplot(plotdf_2, x = 'Genotype', y = 'TPM', color = "black", fill="Genotype",
		title = paste(genetype[i], "pluripotency genes"), add.params = list(size = 1), add =  "dotplot", xlab = " ", palette = EpiLC_XY_palette) +
    	rremove("legend") +
    	stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif") 
	q <- ggpar(q, x.text.angle = 25, ylim = c(0,40), font.main = "bold")
	q <- facet(q, facet.by = "Symbol", nrow = 2)

	TPMplot[[i]] <- q
}

library("gridExtra")
ggsave(snakemake@output[[1]], plot = grid.arrange(grobs = TPMplot, nrow = length(genetype)), width = 8, height = 10)


#plot just 6 of the primed genes for the main figure
primedgenes_6 <-  c("Otx2", "Etv4", "Etv5", "Pou3f1", "Fgf5", "Dnmt3b")

plotdf_6 <- subset(plotdf, type == "primed" & Symbol %in% primedgenes_6)
my_comparisons <- list(c("WT", "5CKO"))
q <- ggboxplot(plotdf_6, x = 'Genotype', y = 'TPM', color = "black", fill="Genotype",
	title = paste("Markers of primed pluripotency"), add =  "dotplot", add.params = list(size = 2), xlab = " ", palette = EpiLC_XY_palette) +
    rremove("legend") +
    stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif") 
q <- ggpar(q, font.main = "bold")
q <- facet(q, facet.by = "Symbol", nrow = 2, scales = "free", panel.labs.font = list(face = 'italic'))


library("gridExtra")
ggsave(snakemake@output[[2]], plot = q, width = 4.5, height = 4.5)


