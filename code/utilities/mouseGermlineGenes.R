#############################################################################################
#                                                                                           #
# 2022 01 24 Determining list of germline-enriched genes for mice                           #
#                                                                                           #
# Inputs:                                                                                   #
#    Mueller_embryonicmouseWWvsequencing.xlsx - file of FPKMs across development of WT and  #
#                 germcell-depleted (WWv) testes                                            #
#    Mueller2013_adultWW_FPKM.txt - file of FPKM from WT WWv.                               #
#                  Raw data from https://pubmed.ncbi.nlm.nih.gov/23872635/                  # 
#    Li_2017_Tissues_FPKMs.xlsx - FPKM across adult mouse tissues                           #
#                  Raw data from https://pubmed.ncbi.nlm.nih.gov/28646208/                  #
#                                                                                           #
#                                                                                           #
# Libraries: readxl, clusterProfiler, org.Mm.eg.db, dyplyr                                  #
# Outputs: Mueller_avgWTgermline_FPKM.txt - table of average fpkm and max fpkm in WT gonads #
#          germGENES20.csv - table of germline-enriched genes with 20% cutoff               #                                                              #
#############################################################################################

#### Goal: find germline-enriched genes by filtering for genes that lose 80% of their expression in the gonads with germ cell depletion

#####1) Find timepoint in WT gonads in which genes of interest have the maximum expression (FPKM)
    #Use this FPKM as the maximum %
    #initially I'm going to use any column with the greatest value, male or female

#####1B) Generate dataframe with WT FPKMs for all gonad sequencing

#file of FPKM values for each embryonic sequencing timepoint
library("readxl")
mueller_genes <- read_excel(snakemake@input[[1]])
mueller_genes <- data.frame(mueller_genes)
#change the tracking ID column to SYMBOL so we can merge together later
names(mueller_genes)[names(mueller_genes) == 'tracking_id'] <- 'SYMBOL'

#File of FPKM values for adult testis
adultfpkm <- read.csv(snakemake@input[[2]], header=TRUE,sep="\t", row.names=1)
#get the average FPKM for the adult timepoints
fpkm_WT = apply(adultfpkm[1:2], 1, mean)
fpkm_WWv = apply(adultfpkm[3:4], 1, mean)
avgFPKM <- cbind(fpkm_WT, fpkm_WWv)
colnames(avgFPKM) <- c("AdultXYWT","AdultXYWWv")
avgFPKM <- as.data.frame(avgFPKM)
avgFPKM <- data.frame(ENSEMBL = rownames(avgFPKM),avgFPKM)

## get gene names with clusterProfiler
library("clusterProfiler")
library(org.Mm.eg.db)
genedf <- bitr(row.names(avgFPKM), fromType = "ENSEMBL", toType = "SYMBOL", 
            OrgDb = org.Mm.eg.db, drop = FALSE)

adultFPKM <- merge(avgFPKM, genedf, by = "ENSEMBL")

# generate a dataframe with WT fpkms across all gonad timepoints (add in adult)
adultFPKM_WT <- subset(adultFPKM, select = c(-AdultXYWWv))
mueller_WT <- subset(mueller_genes, select = c(SYMBOL, P6XYWT, E12XXWT ,E12XYWT, E14XXWT , E14XYWT , E16XXWT , E16XYWT))
mueller_WT <- merge(mueller_WT, adultFPKM_WT, by = "SYMBOL")
#change the column order so the numbers are all together and in chronological order
mueller_WT <- mueller_WT[,c("SYMBOL", "ENSEMBL", "E12XXWT" ,"E12XYWT", "E14XXWT" , "E14XYWT" , "E16XXWT" , "E16XYWT", "P6XYWT", 'AdultXYWT')]
mueller_WT <- unique(mueller_WT)

#print out number of genes with all timepoints that we can use for analysis 
    #19160 genes
cat("Number of genes for analysis: ", nrow(mueller_WT), "\n")

#####1C) Generate dataframe with germ-cell depleted WWv FPKMs for all gonad sequencing
adultFPKM_WWv <- subset(adultFPKM, select = c(-AdultXYWT))
mueller_WWv <- subset(mueller_genes, select = c(SYMBOL, P6XYWWv, E12XXWWv ,E12XYWWv, E14XXWWv , E14XYWWv , E16XXWWv , E16XYWWv))
mueller_WWv <- merge(mueller_WWv, adultFPKM_WWv, by = "SYMBOL")
#change the column order so the numbers are all together and in chronological order
mueller_WWv <- mueller_WWv[,c("SYMBOL", "ENSEMBL", "E12XXWWv" ,"E12XYWWv", "E14XXWWv" , "E14XYWWv" , "E16XXWWv" , "E16XYWWv", "P6XYWWv", "AdultXYWWv")]
mueller_WWv <- unique(mueller_WWv)


#####2) Get the maxiumum WT FPKM value for any timepoint for either sex
mueller_WT[, "WT_max"] <- apply(mueller_WT[, 3:ncol(mueller_WT)], 1, max)


#####3) divide WWv values by WT max value to get the ratio of maximum expression
#make a dataframe that contains the WWv values and the WT max values
WT_max <- subset(mueller_WT, select = c(SYMBOL, ENSEMBL, WT_max))
#keep only genes that have a greater than 1 FPKM in WT_max
WT_max <- subset(WT_max, WT_max > 1)
WWvmaxNORM <- merge(mueller_WWv, WT_max, by = "ENSEMBL")
WWvmaxNORM <- unique(WWvmaxNORM)
WWvmaxNORM <- subset(WWvmaxNORM, select = -c(SYMBOL.y))
#divide each row by WT max
divi <- data.frame(WWvmaxNORM[,1:2], WWvmaxNORM[,3:(ncol(mueller_WT)-1)]/WWvmaxNORM[,ncol(mueller_WT)])
#NaNs are when the value is 0 in WWv. Set NaNs to zero so that we can keep genes that lose expression in WWv.
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
divi[is.nan(divi)] <- 0

#combine together into one df, removing old WT max column
WWvmaxNORM <- merge(WT_max, divi, by = "ENSEMBL")
#drop extra symbol column
WWvmaxNORM <- subset(WWvmaxNORM, select=c(-SYMBOL.x))
print(head(WWvmaxNORM))


#####4) Keep genes that a ratio less than 0.2 (20% of the maximum expression) across all timepoints and sexes
library(dplyr)
muellerRATIO20 = WWvmaxNORM %>% filter_at(vars(4:11), all_vars(. < 0.2))
#this leaves 3800 genes

#####5) Account for expression in other non-gonadal tissues using Li et al 2017 sequencing (https://pubmed.ncbi.nlm.nih.gov/28646208/)
#####5A) Get average FPKMs for Li 2017
Li_genes <- read_excel(snakemake@input[[3]])
Li_genes <- data.frame(Li_genes)

#change into numeric in case imported as character
library(dplyr)
Li_genes_num <- Li_genes[,4:ncol(Li_genes)] %>% mutate_if(is.character,as.numeric)
Li_genes <- data.frame(Li_genes[,1:3], Li_genes_num)

#average expression for males and females in each tissue
Li_Averages <- data.frame(ENSEMBL = Li_genes[,1], SYMBOL = Li_genes[,2], AdrenalGland_F = rowMeans(Li_genes[,4:5]), AdrenalGland_M = rowMeans(Li_genes[,6:7]),
                          Brain_F = rowMeans(Li_genes[,8:11]), Brain_M = rowMeans(Li_genes[,12:15]),
                          Forestomach_F = rowMeans(Li_genes[,16:17]), Forestomach_M = rowMeans(Li_genes[,18:19]),
                          Heart_F = rowMeans(Li_genes[,20:21]), Heart_M = rowMeans(Li_genes[,22:23]),
                          Kidney_F = rowMeans(Li_genes[,24:25]), Kidney_M = rowMeans(Li_genes[,26:27]),
                          Liver_F = rowMeans(Li_genes[,28:31]), Liver_M = rowMeans(Li_genes[,32:35]),
                          
                          LargeIntestine_F = rowMeans(Li_genes[,36:37]), LargeIntestine_M = rowMeans(Li_genes[,38:39]),
                          Lung_F = rowMeans(Li_genes[,40:41]), Lung_M = rowMeans(Li_genes[,42:43]),
                          Muscle_F = rowMeans(Li_genes[,44:45]), Muscle_M = rowMeans(Li_genes[,46:47]),
                          Ovary_F = rowMeans(Li_genes[,48:49]),
                          SmallIntestine_F = rowMeans(Li_genes[,50:51]), SmallIntestine_M = rowMeans(Li_genes[,52:53]),
                          Spleen_F = rowMeans(Li_genes[,54:55]), Spleen_M = rowMeans(Li_genes[,56:57]),
                          Stomach_F = rowMeans(Li_genes[,58:59]), Stomach_M = rowMeans(Li_genes[,60:61]),
                          Testis_M = rowMeans(Li_genes[,62:63]),
                          Thymus_F = rowMeans(Li_genes[,64:65]), Thymus_M = rowMeans(Li_genes[,66:67]),
                          Uterus_F = rowMeans(Li_genes[,68:69]),
                          VesicularGland_M = rowMeans(Li_genes[,70:71]))

#remove testis and ovary since we don't want to include them in the filter
Li_norepro <- subset(Li_Averages, select = -c(Testis_M,Ovary_F))

#ratio of WT testis maximum expression - add the WT max value to the tissues dataset
Li_ratio <- merge(WT_max, Li_norepro, by = "ENSEMBL")
Li_ratio <- subset(Li_ratio, select = -c(SYMBOL.y))
#divide the tissue fpkms by the WT maximum expression
divi_2 <- data.frame(Li_ratio[,1:2], Li_ratio[,4:ncol(Li_ratio)]/Li_ratio$WT_max)
    #NaNs are when the value is 0 in WWv. Set NaNs to zero so that we can keep genes that lose expression in WWv.
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
divi_2[is.nan(divi_2)] <- 0
#add WT max value back
Li_ratio <- merge(WT_max, divi_2, by = "ENSEMBL")
Li_ratio <- subset(Li_ratio, select = c(-SYMBOL.x))

####5B) filter for 20% expression in Li WT tissues 
library(dplyr)
liRATIO20 = Li_ratio %>% filter_at(vars(4:ncol(Li_ratio)), all_vars(. < 0.2))

####6) Subset for genes that are in the Li list AND in the mueller list 
    ### This will give us genes that lose 80% of their maximum expression in WT somatic tissues or with germ cell depletion

germGENES20 <- subset(muellerRATIO20, ENSEMBL %in% liRATIO20$ENSEMBL)
print(head(germGENES20))
germGENES20 <- germGENES20[!(duplicated(germGENES20) | duplicated(germGENES20, fromLast = TRUE)), ]
cat("Number of germline-enriched genes: ", nrow(germGENES20), "\n")
#Number of germline-enriched genes: 2333
germGENES20 <- germGENES20[,1:2]

####7) add in whether the genes are sex biased
# subset muellter dataframe for sex
XXtime <- mueller_WT[,c("SYMBOL", "ENSEMBL", "E12XXWT", "E14XXWT" , "E16XXWT")]
XYtime <- mueller_WT[,c("SYMBOL", "ENSEMBL", "E12XYWT", "E14XYWT" , "E16XYWT")]

Li_ovary <- subset(Li_Averages, select = c(ENSEMBL, Ovary_F))
XXtime <- merge(XXtime, Li_ovary, by = "ENSEMBL")

Li_testis <- subset(Li_Averages, select = c(ENSEMBL, Testis_M))
XYtime <- merge(XYtime, Li_testis, by = "ENSEMBL")

#get the maximum value
XXtime[, "XX_max"] <- apply(XXtime[, 3:ncol(XXtime)], 1, max)
XYtime[, "XY_max"] <- apply(XYtime[, 3:ncol(XYtime)], 1, max)

XXvsXY <- merge(XXtime, XYtime, by = "ENSEMBL")
XXvsXY$XXoverXY <- XXvsXY$XX_max/XXvsXY$XY_max
#using 20% max expression cutoff again
XXvsXY$sexBias <-ifelse(XXvsXY$XXoverXY >= 5,"XX", ifelse(XXvsXY$XXoverXY <= 0.2, "XY", "unbiased"))
XXvsXY_test <- XXvsXY
XXvsXY <- subset(XXvsXY, select = c(ENSEMBL, sexBias))
print(paste0("XXvsXY rows: ", nrow(XXvsXY)))
print(paste0("germ genes not in XX vs XY: ", length(germGENES20$ENSEMBL[!germGENES20$ENSEMBL %in% XXvsXY$ENSEMBL])))


#testing incomplete cases 
germtest <- merge(germGENES20, XXvsXY_test,  by = "ENSEMBL")
incomp <- germtest[!complete.cases(germtest), ]
print("Incomplete cases")
print(nrow(incomp))
print(head(incomp))



germGENES20 <- merge(germGENES20, XXvsXY, by = "ENSEMBL")
germGENES20 <- germGENES20[!(duplicated(germGENES20) | duplicated(germGENES20, fromLast = TRUE)), ]
#keep only complete cases
germGENES20 <- germGENES20[complete.cases(germGENES20), ]

print(paste0("number of germline genes: ", nrow(germGENES20)))


cat("Number of XX germline-enriched genes: ", nrow(subset(germGENES20, sexBias == "XX")), "\n")
cat("Number of XY germline-enriched genes: ", nrow(subset(germGENES20, sexBias == "XY")), "\n")
cat("Number of unbiased germline-enriched genes: ", nrow(subset(germGENES20, sexBias == "unbiased")), "\n")
#XX: 164
#XY: 1195

write.table(germGENES20, file = snakemake@output[[1]], sep =",", row.names = FALSE)

# XXtime[, "XX_avg"] <- apply(XXtime[, 3:ncol(XXtime)], 1, mean)

# XYtime[, "XY_avg"] <- apply(XYtime[, 3:ncol(XYtime)], 1, mean)

# XXvsXY <- merge(XXtime, XYtime, by = "ENSEMBL")
# XXvsXY$XXoverXY <- XXvsXY$XX_avg/XXvsXY$XY_avg
# XXvsXY$sexBias <-ifelse(XXvsXY$XXoverXY > 10,"XX", ifelse(XXvsXY$XXoverXY < 0.1, "XY", "unbiased"))
# checking <- subset(XXvsXY, select = c(ENSEMBL, SYMBOL.x, XX_avg, XY_avg, XXoverXY, sexBias))

# XXvsXY <- subset(XXvsXY, select = c(ENSEMBL, sexBias))


## make a snakey graph of filtering for germline-enriched genes
#https://github.com/davidsjoberg/ggsankey
#devtools::install_github("davidsjoberg/ggsankey")