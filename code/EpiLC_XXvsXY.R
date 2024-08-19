### 23.09.08 euler diagram comparing XX and XY EpiLCs with loss of Kdm5c

##1) Read in the germline genes list
germ <- read.csv(snakemake@input[[1]], sep = ",")

source("code/utilities/parameters.R") #get the log2fc cut off

##2) Read in the EpiLC XX and XY germline DEGs
samples <- c("XY5cKO", "XX5cHET", "XX5cKO")
DEGs <- list()
for (s in 1:length(samples)){
	#read in the dataframe
	df <- read.csv(snakemake@input[[s+1]], sep = ",")
	
	gdf <- subset(df, ENSEMBL %in% germ$ENSEMBL & log2FoldChange > l2fcco_ESCEpiLC)
	germGENES <- gdf$ENSEMBL
	#subset for the number of germline genes
	DEGs[[s]] <- germGENES
}
names(DEGs) <- samples


##3) Calculate the overlap between groups
#overlap <- sapply(DEGs, function(x) sapply(DEGs, function(y) sum(y %in% x)))
#print(overlap)

#genes in all groups
all <- Reduce(intersect, DEGs)

#genes in pair combinations
paired_overlap <- function(a, b){
	ol <- intersect(a, b)
	ol_unique <- ol[!ol %in% all]
	return(ol_unique)
}

XY5cKO_AND_XX5cHET <- paired_overlap(DEGs[["XY5cKO"]], DEGs[["XX5cHET"]])
XY5cKO_AND_XX5cKO <- paired_overlap(DEGs[["XY5cKO"]], DEGs[["XX5cKO"]])
XX5cKO_AND_XX5cHET <- paired_overlap(DEGs[["XX5cKO"]], DEGs[["XX5cHET"]])

#genes only in each sample
#a is the one you're keeping, b and c are the ones you're comparing to
unique_overlap <- function(a, b, c){
	unique <- a[!a %in% b & !a %in% c]
	return(unique)
}

XY5cKO_only <- unique_overlap(DEGs[["XY5cKO"]], DEGs[["XX5cHET"]], DEGs[["XX5cKO"]])
XX5cHET_only <- unique_overlap(DEGs[["XX5cHET"]], DEGs[["XY5cKO"]], DEGs[["XX5cKO"]])
XX5cKO_only <- unique_overlap(DEGs[["XX5cKO"]], DEGs[["XX5cHET"]], DEGs[["XY5cKO"]])


##4) Make a dataframe with each gene in the list and which is female or male specific

EpiLC_germ_XXvsXY <- data.frame(ENSEMBL = unlist(DEGs))
print(head(EpiLC_germ_XXvsXY))

#which sample it is a DEG
EpiLC_germ_XXvsXY$Sample_DEG <- ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XY5cKO_only, "XY 5CKO",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XY5cKO_only, "XY 5CKO", 
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XX5cHET_only, "XX 5CHET",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XX5cKO_only, "XX 5CKO",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  all, "All",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XY5cKO_AND_XX5cHET, "XY 5CKO and XX 5CHET",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XY5cKO_AND_XX5cKO, "XY 5CKO and XX 5CKO",
	ifelse(EpiLC_germ_XXvsXY$ENSEMBL %in%  XX5cKO_AND_XX5cHET, "XX 5CKO and XX 5CHET", "whoops"
	))))))))

#log2fold change of the gene in each genotype

EpiLC_germ_XXvsXY <- unique(merge(EpiLC_germ_XXvsXY, germ, by = "ENSEMBL"))

print("EpiLC_germ_XXvsXY")
print(head(EpiLC_germ_XXvsXY))
write.table(EpiLC_germ_XXvsXY, snakemake@output[[1]], sep = ",", row.names = F)



##4) Make the euler plot
library("eulerr")
source("code/utilities/colorpalettes.R")
together <- euler(c("XY5cKO" = length(XY5cKO_only), "XX5cHET" = length(XX5cHET_only), "XX5cKO" = length(XX5cKO_only),
                "XY5cKO&XX5cHET" = length(XY5cKO_AND_XX5cHET), "XY5cKO&XX5cKO" = length(XY5cKO_AND_XX5cKO), "XX5cHET&XX5cKO" = length(XX5cKO_AND_XX5cHET),
                "XY5cKO&XX5cHET&XX5cKO" = length(all)))

p <- plot(together, quantities = TRUE, labels = list(font = 4), fills = c(EpiLC_XY_KO, EpiLC_XX_HET, EpiLC_XX_KO))

library("ggplot2")
ggsave(snakemake@output[[2]], plot = p, width = 4, height = 4)


#24.03.04 Make a histogram of the number of unique vs shared germline genes per chromosome to show XX increase is not an X-inactivation defect
#also make a graph of just the total number per chromosome?
#should it be absolute number of % of all germline genes on that chromosome?

#1) get the chromosome location of the all germline genes
source("code/utilities/germ_chromo.R")

#2) get which germline genes are unique between males and females
XXallDEGs <- c(DEGs[["XX5cHET"]], DEGs[["XX5cKO"]])
XX_UNIQUE <- setdiff(XXallDEGs, DEGs[["XY5cKO"]])
XY_UNIQUE <- setdiff(DEGs[["XY5cKO"]], XXallDEGs)

#they are in list format
xxonly_chr <- chromo_DEG_count_list(XX_UNIQUE)
shared_chr <- chromo_DEG_count_list(all)

print("xx only:")
print(xxonly_chr)
print("shared:")
print(shared_chr)

#plot the percentage (germ genes/ all chromosome genes)
	#chr_cnt_df - dataframe of the chromo counts (outdf from above)
	#y_col_name - name of the column you want for the y-axis, as string, that way same code can be used for multiple
	#title_name - title of plot
	#out - number of snakemake output for saving plot
chromo_hist(shared_chr, "DEG_ratio", "shared germline DEGs", 3)
chromo_hist(xxonly_chr, "DEG_ratio", "XX only germline DEGs", 4)

#save the chromosome results
write.table(xxonly_chr, snakemake@output[[5]], sep = ",", row.names = F)

### histogram of total number of germline genes
source('code/utilities/colorpalettes.R')
print("checking samples:")
print(samples)
count <- c()
for (i in 1:length(samples)){
	count[i] <- length(DEGs[[i]])
}

germhist <- data.frame(samples, count)
germhist$samples <- factor(germhist$samples, levels = samples)
print(germhist)

library('ggpubr')
q <- ggbarplot(germhist, x = "samples", y = "count", color = "samples", fill = "samples", palette = EpiLCpalette, label = TRUE, lab.pos = "in", lab.col = "black")

ggsave(snakemake@output[[6]], plot = q, width = 3, height = 3)


### Gene ontology of clusters
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)

sharedgroups <- c("All", "XY 5CKO and XX 5CHET", "XY 5CKO and XX 5CKO")
EpiLC_germ_XXvsXY$category <- ifelse(EpiLC_germ_XXvsXY$Sample_DEG %in% sharedgroups, "Shared", EpiLC_germ_XXvsXY$Sample_DEG)

EpiLC_germ_XXvsXY$category[EpiLC_germ_XXvsXY$category == "XY 5CKO"] <- "Male DEG"
EpiLC_germ_XXvsXY$category[EpiLC_germ_XXvsXY$category == "XX 5CHET" | EpiLC_germ_XXvsXY$category == "XX 5CKO" | EpiLC_germ_XXvsXY$category == "XX 5CKO and XX 5CHET"] <- "Female DEG"

#keep only the ensembl and cluster info
clusterDEGs <- subset(EpiLC_germ_XXvsXY, select = c(ENSEMBL, category))
#print(head(clusterDEGs))

#split the dataframe into a list of dataframes based on cluster
germDEGs <- list()
#for every cluster, get the ENSEMBL
for (i in unique(clusterDEGs$category)){
	sub <- subset(clusterDEGs, category == i)
	germDEGs[[i]] <- sub$ENSEMBL
}

#germDEGs <- split(clusterDEGs , f = clusterDEGs$Cluster)
names(germDEGs) <- unique(clusterDEGs$category)

#now run the gene ontology comparison
ck <- compareCluster(geneCluster = germDEGs, fun = enrichGO,  OrgDb = "org.Mm.eg.db", keyType="ENSEMBL", ont="BP")
#ck <- setReadable(ck, OrgDb = "org.Mm.eg.db", keyType="ENSEMBL")
head(ck) 
ck <- simplify(ck)

write.table(ck, snakemake@output[[7]], row.names = FALSE, sep = ",")

ggsave(snakemake@output[[8]], plot = dotplot(ck), width = 5.5, height = 5.5)





