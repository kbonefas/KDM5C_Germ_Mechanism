#############################################################################################
#                                                                                           #
# 2022 11 08 For pulling out DEGs from a DESeq2 results table using Snakemake                              # 
#                                                                                           #
# Inputs: DESeq2 results file (.csv), padj cut off                                          #
# Libraries: clusterProfiler, org.Mm.eg.db                                                            #
# Outputs: csv file of upregulated and downregulated DEGs                                   #
#                                                                                           #
#############################################################################################
###function to pull out DEGs
    #alph - alpha cut off value for significance for DESeq2
    #n = iteration for for loop

source("code/utilities/parameters.R") #get cut off of log2fc (l2fcco)

DEGtable <- function(alph, n){
    res <- read.csv(snakemake@input[[n]], sep =",", row.names = 1)
        #make a new column for ensembl names
    resdf <- data.frame(res, ENSEMBL = row.names(res))
        #sort out the DEGs
    resDEG <- subset(resdf, padj < alph & abs(log2FoldChange) > abs(l2fcco))
        #add a column for the direction of the DEG (Downregulated or Upregulated)
    resDEG$Direction <- ifelse(resDEG$log2FoldChange > l2fcco, "Up", "Down")

        #use clusterProfiler to covert ensembl IDs to gene symbols 
    library("clusterProfiler")
    library(org.Mm.eg.db)
    genedf <- bitr(row.names(resDEG), fromType = "ENSEMBL", toType = "SYMBOL", 
                        OrgDb = org.Mm.eg.db, drop = FALSE) #if the ids do not have a gene symbol, do not drop them from the table

        #merge gene name conversion with DEG dataframe
    resDEG_ids <- merge(resDEG, genedf, by = "ENSEMBL")
	#order by log2fc and padj
	resDEG_ids_sort <- resDEG_ids[order(-resDEG_ids$log2FoldChange, resDEG_ids$padj),]
    print(head(resDEG_ids_sort))
        #write the table
    write.table(resDEG_ids_sort, snakemake@output[[n]], sep=",", row.names = FALSE)
}

    #### loop for each comparison
# for (i in 1:snakemake@params[["numbcompare"]]){
#     DEGtable(snakemake@params[["alpha"]], i)
# }

#for use when res table is already in R
#res - results table
#n - position for output
DEGtable2 <- function(res, n, alph){
        #make a new column for ensembl names
    resdf <- data.frame(res, ENSEMBL = row.names(res))
        #sort out the DEGs
    resDEG <- subset(resdf, padj < alph & abs(log2FoldChange) > abs(l2fcco))
        #add a column for the direction of the DEG (Downregulated or Upregulated)4
    resDEG$Direction <- ifelse(resDEG$log2FoldChange > l2fcco, "Up", "Down")

        #use clusterProfiler to covert ensembl IDs to gene symbols 
    library("clusterProfiler")
    library(org.Mm.eg.db)
    genedf <- bitr(row.names(resDEG), fromType = "ENSEMBL", toType = "SYMBOL", 
                        OrgDb = org.Mm.eg.db, drop = FALSE)#if the ids do not have a gene symbol, do not drop them from the table

        #merge gene name conversion with DEG dataframe
    resDEG_ids <- merge(resDEG, genedf, by = "ENSEMBL")
	#order by log2fc and padj
	resDEG_ids_sort <- resDEG_ids[order(-resDEG_ids$log2FoldChange, resDEG_ids$padj),]
    #print(head(resDEG_ids_sort))
        #write the table
    write.table(resDEG_ids_sort, snakemake@output[[n]], sep=",", row.names = FALSE)
	return(resDEG_ids_sort)
}