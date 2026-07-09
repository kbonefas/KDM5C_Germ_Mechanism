#2026.07.09 - chromosome distribution of germline genes 
    # are certain types of germline genes enriched on certain chromosomes?


#sperm, egg, unbiased
#CGI vs non-CGI
#is there overall enrichment on the X chromosome - https://www.nature.com/articles/ng.2705, https://www.nature.com/articles/ng.126


#1) read in csv file with the gene names
germ <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)

source("code/utilities/GeneTSSandTES.R")
library("GenomicRanges")

GR_list <- list()

#for loop to iterate through the categories
print(unique(germ$sexBias))
for(i in 1:length(unique(germ$sexBias))){
    print(unique(germ$sexBias)[i])
    # subset the germ dataframe for the current category
    germ_subset <- subset(germ, germ$sexBias == unique(germ$sexBias)[i])

    #2) get the gene coordinates
    gene_coords <- geneTSSandTES_df(goi = germ_subset$ENSEMBL)

    #3) make a genomic ranges object for the genes
    query <- makeGRangesFromDataFrame(gene_coords, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "TSS", end.field = "TES")
	print(head(query))

    #4) save the genomic ranges object to a list
    GR_list[[i]] <- query
}

#6) read into genomicDistributions


library("GenomicDistributions")

# First, calculate the distribution:
queryList <- GRangesList(unbiased=GR_list[[3]], sperm=GR_list[[1]], egg=GR_list[[2]])
x2 <- calcChromBinsRef(queryList, "mm10")

#7) plot

pdf(file = snakemake@output[[1]], width = 10, height = 10)
	plotChromBins(x2)
dev.off()


