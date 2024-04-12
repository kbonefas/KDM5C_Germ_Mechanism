####################################
# 11.28.2022 
# MA plot of imprinted genes 
#
#####################################
library(ggpubr)
#read in the imprinted gene list from https://www.geneimprint.com/site/genes-by-species.Mus+musculus
impgenes <- read.csv("data/raw/220610_geneimprint_mm.csv", sep=',')
#keep only the imprinted genes that are cannonically imprinted and expressed from either maternal or paternal alleles
impgenes <- subset(impgenes, impgenes$Status == "Imprinted")
impgenes <- subset(impgenes, impgenes$Expressed.Allele == "Paternal" | impgenes$Expressed.Allele == "Maternal")


impMA <- function(res){
	#make new column called gene with the ensmbl ids
	res$Gene <- row.names(res)
	res$log2baseMean <- log2(res$baseMean)
	#merge imprinted list with res table so we have expression of imprinted genes
	impres <- merge(res, impgenes, by = "Gene")
	#get which genes are the degs so we can label only those
	degs <- subset(impres, impres$padj < snakemake@params[["alpha"]])

	p <- ggscatter(impres, x = "log2baseMean", y = "log2FoldChange", color = "Expressed.Allele", palette = c("#db72a3", "#639fbf"),
		label = "Symbol", repel = TRUE, label.select = list(criteria = "'padj' < 0.1"))
	ggpar(p, ylim = c(-7,7))
} 

#read in the DESeq2 res table
for (i in 1:snakemake@params[["numbcompare"]]){
	res <- read.csv(snakemake@input[[i]], sep =",", row.names = 1)
	p <- impMA(res)
	ggsave(snakemake@output[[i]], plot = p, width = 10, height = 10)

}





