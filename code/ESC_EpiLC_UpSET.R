#2023.01.13 - Upset plot germline genes 

#read in the DEGs
samples <- c("ESC", "48VA", "48NO", "96VA", "96NO")
for(i in 1:length(samples)){
    DEGs <- read.csv(snakemake@input[[1+i]], sep = ",", row.names = 1)
    print(head(DEGs))
}

germ <- read.csv(snakemake@input[[1]])
print("germline genes")
print(head(germ))