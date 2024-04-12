#############################################################################################
#                                                                                           #
# 2022 10 25 Enrichplot for gene ontology analysis                                          #
#                                                                                           #
# Inputs: DEGs results table (.csv)                                                         #
# Libraries: enrichplot, org.Mm.eg.db, clusterProfiler, ggplot2                             #
# Outputs: Plot of gene ontologies in a network map                                         #
# Citation: Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler:           #
#   an R package for comparing biological themes among gene clusters.                       #
#   OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287                            #
#                                                                                           #
#############################################################################################

#degs - dataframe of DEGs (generated after running DEGtable2 from code/utilities/DESeq2_DEGs.R)
#FCcutoff - what is the cut off value for the fold change?
GOmap <- function(degs, FCcutoff){
    library(enrichplot)
    library(org.Mm.eg.db)
    library(clusterProfiler)
    library(ggplot2)


	#make a dataframe with the first column ENSEMBL IDs and the second column log2FoldChange
		#first and third columns of the degs table
	#df <- data.frame(ENSEMBL = degs[,1], log2FoldChange = degs[,3])

	#sort the genes log2fc from low to high and set the names to the gene IDs
    geneList <- degs[,3]
    names(geneList) <- as.character(degs[,1])
    geneList <- sort(geneList, decreasing = TRUE)

	#get the names that meet the the cut-off value for log2fc
    de <- names(geneList)[geneList > FCcutoff]
    #choose what organism and what type of ontology (biological process)
    ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
    
	edo <- pairwise_termsim(ego)
	p <- emapplot(edo)
	return(p)

    #need to make 1st column Entrez gene ID and 2nd column fold change
    # #convert to entrez
    # entrez <- bitr(degs$ENSEMBL, fromType = "ENSEMBL",
    #                 toType = c("ENTREZID"),
    #                 OrgDb = org.Mm.eg.db, drop = FALSE)
    #merge back with the original df and subset for needed columns
    # entrez <- merge(entrez, DEGs, by = "ENSEMBL")                
    # entrez <- subset(entrez, select =c(ENTREZID, log2FoldChange))

    #sort the log2fc based on decreasing to increasing and set the names to the gene IDs
    # geneList <- entrez[,2]
    # names(geneList) <- as.character(entrez[,1])
    # geneList <- sort(geneList, decreasing = TRUE)

    #choose cut-off value for log2fc
    # de <- names(geneList)[geneList > FCcutoff]
    # #choose what organism and what type of ontology (biological process)
    # ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
    
    # ggsave(paste0("results/GOmap_", title, "_fc", FCcutoff, ".pdf"), plot = pl, width = 7, height = 7)


}

GOdot <- function(degs, FCcutoff){
    library(enrichplot)
    library(org.Mm.eg.db)
    library(clusterProfiler)
    library(ggplot2)


	#make a dataframe with the first column ENSEMBL IDs and the second column log2FoldChange
		#first and third columns of the degs table
	#df <- data.frame(ENSEMBL = degs[,1], log2FoldChange = degs[,3])

	#sort the genes log2fc from low to high and set the names to the gene IDs
    geneList <- degs[,3]
    names(geneList) <- as.character(degs[,1])
    geneList <- sort(geneList, decreasing = TRUE)

	#get the names that meet the the cut-off value for log2fc
    de <- names(geneList)[geneList > FCcutoff]
    #choose what organism and what type of ontology (biological process)
    ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
	print(ego[, c("ID", "Description", "p.adjust")])

	
	q <- dotplot(ego, showCategory=15)
	return(q)

    #need to make 1st column Entrez gene ID and 2nd column fold change
    # #convert to entrez
    # entrez <- bitr(degs$ENSEMBL, fromType = "ENSEMBL",
    #                 toType = c("ENTREZID"),
    #                 OrgDb = org.Mm.eg.db, drop = FALSE)
    #merge back with the original df and subset for needed columns
    # entrez <- merge(entrez, DEGs, by = "ENSEMBL")                
    # entrez <- subset(entrez, select =c(ENTREZID, log2FoldChange))

    #sort the log2fc based on decreasing to increasing and set the names to the gene IDs
    # geneList <- entrez[,2]
    # names(geneList) <- as.character(entrez[,1])
    # geneList <- sort(geneList, decreasing = TRUE)

    #choose cut-off value for log2fc
    # de <- names(geneList)[geneList > FCcutoff]
    # #choose what organism and what type of ontology (biological process)
    # ego <- enrichGO(de, keyType = 'ENSEMBL', OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
    
    # ggsave(paste0("results/GOmap_", title, "_fc", FCcutoff, ".pdf"), plot = pl, width = 7, height = 7)


}
