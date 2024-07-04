
## 2024.06.30 Get the transcript start and end site of each gene in the clusters for bedtools plotting of KDM5C binding
source("code/utilities/GeneTSSandTES.R")

#tss window
WIND <- 1000

#germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")
geneTSSandTES(germ$ENSEMBL, 1)
geneTSSwindow(germ$ENSEMBL, 2, WIND)

#test genes 
testgenes <- c("ENSMUSG00000010592", "ENSMUSG00000020059", "ENSMUSG00000009628")
geneTSSandTES(testgenes, 3)
geneTSSwindow(testgenes, 4, WIND)


#just the brain or EpiLC DEGs
samples <- c("EpiLC", "Amygdala", "Hippocampus")
allgermDEGs <- c()
allupDEGs <- c()
for (i in 1:length(samples)){
	DEGs <- read.csv(snakemake@input[[i+1]], sep = ",")
	#non-germline up DEGs
	upDEGs <-subset(DEGs, Direction == "Up" & !(ENSEMBL %in% germ$ENSEMBL))
	allupDEGs <- c(allupDEGs, upDEGs$ENSEMBL)

	#germline up DEGs
	germDEGs <-subset(DEGs, Direction == "Up" & ENSEMBL %in% germ$ENSEMBL)
	allgermDEGs <- c(allgermDEGs, germDEGs$ENSEMBL)

}

allgermDEGs <- unique(allgermDEGs)
print(length(allgermDEGs))
geneTSSandTES(allgermDEGs, 5)
geneTSSwindow(allgermDEGs, 6, WIND)


allupDEGs <- unique(allupDEGs)
print(length(allupDEGs))
geneTSSandTES(allupDEGs, 7)
geneTSSwindow(allupDEGs, 8, WIND)
