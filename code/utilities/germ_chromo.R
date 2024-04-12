### 23.12.12 - Visualize chromosome location of germline DEGs
print(getwd())
#1) Read in the germline genes

germ <- read.csv("data/processed/germGENES20.csv")

#2) get the chromosome location of the all germline genes

#read in gtf file of all mm10 transcripts
#Using same gtf as the alignment, mm10 (GRCm38.p6) from gencode - https://www.gencodegenes.org/mouse/release_M23.html

library("rtracklayer")
gtf <- readGFF("data/raw/gencode.vM23.annotation.gtf")
gtf_df <- as.data.frame(gtf)
print(head(gtf_df))

#get just the geneID, split transcript variant into new column so you just have ensembl name
ENSEMBL_chr <- subset(gtf_df, select = c("gene_id", "seqid"))
ENSEMBL<- read.csv(text = ENSEMBL_chr$gene_id, sep=".", header = FALSE, col.names = c("gene_id", "var"))
ENSEMBL$chr <- ENSEMBL_chr$seqid
names(ENSEMBL)[names(ENSEMBL) == 'gene_id'] <- 'ENSEMBL'
#remove var column
ENSEMBL2 <- ENSEMBL[,-2]
#keep only the unique columns (so there is only one entry per gene)
ENSEMBL2 <- unique(ENSEMBL2)

#merge with germline gene list to get chromosome locations
germ_chromo <- merge(germ, ENSEMBL2, by= 'ENSEMBL')
print(tail(germ_chromo))
print(paste0("germ rows:", nrow(germ), " new rows:", nrow(germ_chromo)))



#3) get the % of germline genes on each chromosome 
#list of chromosomes we're looking at
numbies <- c(1:22) #all of the numbers
print("chromlist:")
chromlist <- unique(ENSEMBL2$chr) #get the unique ones
#order(chromlist, decreasing = FALSE)
print(chromlist)
germ_chromo$chr <- as.character(germ_chromo$chr)



#count the number of all ENSEMBL transcripts on that chromosome
#make a dataframe with each chromosome, the total number of genes, the number of germline genes, and the ratio
total_cnt <- c()
germ_cnt <- c()

for (i in 1:length(chromlist)){
	testchr <- chromlist[i]
	total_cnt[i] <- nrow(subset(ENSEMBL2, chr == testchr))
	germ_cnt[i] <- nrow(subset(germ_chromo, chr == testchr))
}

chromodf <- data.frame(chr = chromlist, total = total_cnt, germline = germ_cnt, ratio = germ_cnt/total_cnt)
print(chromodf)


#4) calculate number and % of germline DEGs on chromosome
source("code/utilities/parameters.R") #get the log2fc cut off
library("ggplot2")

#DEGs - dataframe of all DEGs for a sample, with ENSEMBL IDs and log2FoldChange
chromo_DEG_count <- function(DEGs){
	##1) read in the DEGs
	DEGs_up <- subset(DEGs, log2FoldChange > l2fcco)

	##2) subset for germline genes by merging upregulated genes with germ list
	DEGs_germ_chromo <- merge(DEGs, germ_chromo, by = "ENSEMBL")

	##3) get the # of germline DEGs on each chromosome 
	DEG_cnt <- c()
	for (i in 1:length(chromlist)){
	testchr <- chromlist[i]
	DEG_cnt[i] <- nrow(subset(DEGs_germ_chromo, chr == testchr))}

	#output dataframe with number of germline DEGs on each chromosome and % of all germline genes on the chromosome
	outdf <- data.frame(chromodf, DEGs = DEG_cnt, DEG_ratio = DEG_cnt/germ_cnt)
	return(outdf)
}

#same as above but for list of ensembl genes
chromo_DEG_count_list <- function(genies){
	##1) read in the DEGs

	##2) subset for germline genes by merging upregulated genes with germ list
	DEGs_germ_chromo <- subset(germ_chromo, ENSEMBL %in% genies)

	##3) get the # of germline DEGs on each chromosome 
	DEG_cnt <- c()
	for (i in 1:length(chromlist)){
	testchr <- chromlist[i]
	DEG_cnt[i] <- nrow(subset(DEGs_germ_chromo, chr == testchr))}

	#output dataframe with number of germline DEGs on each chromosome and % of all germline genes on the chromosome
	outdf <- data.frame(chromodf, DEGs = DEG_cnt, DEG_ratio = DEG_cnt/germ_cnt)
	return(outdf)
}

#plot the germline genes
#chr_cnt_df - dataframe of the chromo counts (outdf from above)
#y_col_name - name of the column you want for the y-axis, as string, that way same code can be used for multiple
#title_name - title of plot
#out - number of snakemake output for saving plot
library('ggpubr')
chromo_hist <- function(chr_cnt_df, y_col_name, title_name, out){

	p <- ggbarplot(chr_cnt_df, "chr", y_col_name, fill = "chr", color = "chr") + rotate_x_text(45)
	pl <- ggpar(p, legend = "none", ylab = "# of genes/total", xlab = "Chromosome #", title = paste0("Location of ", title_name) )
	ggsave(snakemake@output[[out]], plot = pl, width = 5, height = 4)

}
### Instead of doing each individually should just add chromosome to the germline genes list info