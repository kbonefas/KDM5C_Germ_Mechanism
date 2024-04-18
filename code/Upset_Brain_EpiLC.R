### 24.03.06 adapted from dissertation code 23.10.20 upset plot comparing XX and XY EpiLCs and brain germline DEGs

##1) Read in the germline genes list
germ <- read.csv(snakemake@input[[1]], sep = ",")

##2) Read in the EpiLC and brain XX and XY germline DEGs, make sure files are in same order
source("code/utilities/parameters.R") #for the log2fc cut off (l2fcco)
samples <- c("EpiLC XY 5CKO", "AMY 5CKO", "HIP 5CKO")
DEGs <- list()
for (s in 1:length(samples)){
	#read in the dataframe
	df <- read.csv(snakemake@input[[s+1]], sep = ",")
	
	gdf <- subset(df, ENSEMBL %in% germ$ENSEMBL & log2FoldChange > l2fcco)
	germGENES <- gdf$ENSEMBL
	#subset for the number of germline genes
	DEGs[[s]] <- germGENES
}
names(DEGs) <- samples

##3) Calculate the overlap between groups

#make an upset plot for the overlap
library("UpSetR")

modifiedupset <- function(samplelist){
	upset(fromList(samplelist), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "# of Overlapping Germline DEGs", empty.intersections = "on")
}



#just male samples
pdf(file = snakemake@output[[1]], width = 8, height = 5)

upset(fromList(DEGs), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "Overlapping Germline DEGs", empty.intersections = "on", text.scale = 2)

dev.off()


