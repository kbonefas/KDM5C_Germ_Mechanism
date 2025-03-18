#23.10.10 Visualizing germline genes dysregulated with Kdm5c loss in ESCs and EpiLCs

### Gene ontology of clusters
clustergenes <- read.csv(snakemake@input[[1]], sep = ",")
print(head(clustergenes))
#24.06.07 gene ontology of 5CKO germline genes in different clusters
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)


clusterDEGs <- clustergenes[,1:2]
#print(head(clusterDEGs))

#split the dataframe into a list of dataframes based on cluster
cluster_GO <- list()
#for every cluster, get the ENSEMBL
for (i in unique(clusterDEGs$Cluster)){
	sub <- subset(clusterDEGs, Cluster == i)
	cluster_GO[[i]] <- sub$ENSEMBL
}

#germDEGs <- split(clusterDEGs , f = clusterDEGs$Cluster)
names(cluster_GO) <- c("1", "2", "3", "4", "5", "6", "7", "8")
print(head(cluster_GO))


#now run the gene ontology comparison
ck <- compareCluster(geneCluster = cluster_GO, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
ck <- simplify(ck) 

write.table(ck, snakemake@output[[1]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[2]], plot = dotplot(ck), width = 7, height = 8.5)


