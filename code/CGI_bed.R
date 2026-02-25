#2026.02.24 - bed files of germline genes with and without CGIs

source("code/utilities/GeneTSSandTES.R")

#tss window
WIND <- 3000

#germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")

germ_CGI <- subset(germ, germ$Promo_CGI == "CGI")
germ_CGIfree <-  subset(germ, germ$Promo_CGI == "no")

geneTSSwindow(germ_CGI$ENSEMBL, 1, WIND)
geneTSSwindow(germ_CGIfree$ENSEMBL, 2, WIND)