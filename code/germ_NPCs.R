#germ germline genes from npcs
source("code/utilities/DESeq2_DEGs.R")

DEGtable(snakemake@params[["alpha"]], 1)