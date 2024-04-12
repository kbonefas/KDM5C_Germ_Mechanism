# 23 01 10 - Get list of DEGs that are germline genes 

#read in the generated germline genes
germ <- read.csv(snakemake@input[[1]], sep =",")

###function to pull out germline DEGs from DEG list
    #n = iteration of for loop (which sample)
germDEGtable <- function(n){
	#start with n+1 because the first input is the germline gene list
    DEGs <- read.csv(snakemake@input[[n+1]], sep =",")
        #sort out the germline DEGs
	germDEGs <- subset(DEGs, DEGs$Direction == "Up")
    germDEGs <- subset(germDEGs, germDEGs$ENSEMBL %in% germ$ENSEMBL)

    write.table(germDEGs, snakemake@output[[n]], sep=",", row.names = FALSE)
}

#function to return dataframe of germline DEGs
#DEGs <- dataframe of DEGs

germDEGs <- function(DEGs){
    #germ <- dataframe of germline genes
    germ <- read.csv("data/processed/germGENES20.csv", sep =",")
    #sort out the germline DEGs
	germDEGs <- subset(DEGs, DEGs$Direction == "Up")
    germDEGs <- subset(germDEGs, germDEGs$ENSEMBL %in% germ$ENSEMBL)

    return(germDEGs)

}


#### loop for each sample
for (i in 1:snakemake@params[["samplenumber"]]){
    germDEGtable(i)
}
