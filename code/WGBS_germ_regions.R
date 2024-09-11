### 24.08.11 - Getting regions downstream from the TSS to compare for WGBS
source("code/utilities/GeneTSSandTES.R")

#tss window - how many basepairs up and downstream downstream from the TSS
WIND <- 500

#germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")
geneTSSwindow(germ$ENSEMBL, 1, WIND)

#all gene promoters using EpiLC gene expected counts Ensembl list
all_genes <- read.csv(snakemake@input[[2]], sep = ",")
geneTSSwindow(all_genes$gene_id, 1, WIND)
