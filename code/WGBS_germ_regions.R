### 24.08.11 - Getting regions downstream from the TSS to compare for WGBS
source("code/utilities/GeneTSSandTES.R")

#tss window - how many basepairs up and downstream downstream from the TSS
WIND <- 500

#germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")
geneTSSwindow(germ$ENSEMBL, 1, WIND)

