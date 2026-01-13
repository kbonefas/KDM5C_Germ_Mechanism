#2023.01.13 - Upset plot germline genes 

#read in the DEGs
samples <- c("ESC", "48VA", "48NO", "96VA", "96NO")
for(i in samples){
    DEGs <- read.csv(snakemake@input[[1+i]], sep = ",", rownames = TRUE)
    print(head(DEGs))
}