#############################################################################################
#                                                                                           #
# 2022 11 7 DESeq2 using snakemake                                                          #
#  This script assumes the only variable you're comparing between samples is their genotype #                                                                                         #
# it also assumes your wildtype sample is called WT
# Inputs:
#       cts = counts file (.txt)                                                            #
#       sampleinfo = sample information sheet #MAKE SURE THE SAMPLES MATCH THE CTS FILE NAMES
# Libraries: DESeq2, ggplot2                                                                #
# Outputs: PCA plot of samples and DESeq2 results tables (.csv)                             #
#                                                                                           #
#############################################################################################

###1) format the raw counts table so that all values are integers and there are no gene variant numbers
cts <- read.csv(snakemake@input[[1]], sep ="\t", row.names = 1)
    #convert counts to integers
cts[1:ncol(cts)] <- lapply(cts[1:ncol(cts)], as.integer)


###2) read in the sample information
    #Note: your sample names MUST match the counts file column names
    #Note: your wildtype samples must have a genotype of WT 
SampleInfo <- read.csv(snakemake@input[[2]], sep =",") 
colnames(SampleInfo) <- c("Sample", "Genotype", "Region", "Type")

###3) Set up and run DESeq2 
library('DESeq2')
    #create the coldata
coldata <- data.frame(row.names = SampleInfo$Sample, genotype = SampleInfo$Genotype, region = SampleInfo$Region, type = SampleInfo$Type)
    #compare all samples
SAMPLES <- coldata$type == "paired.end"
countTable <- cts[ , SAMPLES]
    #Assign the genotype to each column 
genotype <- coldata$genotype[ SAMPLES]

    #make a DESeqDataSet (dds) from the counts table, comparing between genotypes (design)
dds <- DESeqDataSetFromMatrix(countData = countTable,
                            colData = coldata,
                            design = ~ genotype)

    #set WT as the default level for comparison 
dds$genotype <- relevel(dds$genotype, ref = "WT")
    #run DESeq2

    #add in region so that we can account for the interaction between genotype and brain region 
dds$group <- factor(paste0(dds$region, dds$genotype))
    #new design using the group variable
#dds$group <- relevel(dds$group, ref = "amygdalaWT")
design(dds) <- ~ group

    #run deseq on the DESeqDataSet to generate all the results tables
dds <- DESeq(dds)
        #This prints out the coef contrasts created by DESeq2 with WT as the comparison
    #They will be in numerical then alphabetical order for each of your genotypes 
        #Example: 5cKO
        #MAKE SURE the snakefile output data tables matches this order
resultsNames(dds)

###4) generate the PCA plot comparing between groups 
library('ggplot2')
    #variance stabilizing transformation
        #essentially accounts for sampling variablity of low counts while also normalizing variability across all samples 
vsd <- vst(dds, blind=TRUE)
data <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
    #adding region as a separate variable back in so that the shape of the point is dependent upon sex
data <- data.frame(data, Region = SampleInfo$Region)
q <- ggplot(data, aes(PC1, PC2, color=group, shape = Region)) +
    geom_point(size=4.5) +
    ggtitle("PCA Plot of Adult Brain Samples")+
    scale_color_manual(values=c("#E36752", "#732113", "#556FAF","#2D3C61")) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) 

    #Save the PCA plot
ggsave(snakemake@output[[1]], plot = q, width = 5, height = 4)


###5) generate results tables for each comparision  
    #This function runs the results and the shrink for each comparison
DESeq2_restab <- function(n, i){
    res_unshrunk <- results(dds, contrast = n, alpha = snakemake@params[["alpha"]])
    #LFCshrink to normalize small count variation
    res <- lfcShrink(dds, contrast=n, res=res_unshrunk, type = "ashr")
    #prints the results summary (# of DEGs)
    summary(res)
    write.table(res, snakemake@output[[i]], sep=",", col.names=NA)
    return(res)
}
 

AMY_5CWT <- DESeq2_restab(c("group", "amygdala5cKO", "amygdalaWT")   , 2)
#change the reference level to the wild type sample that you need
# dds$group <- relevel(dds$group, ref = "hippocampusWT")
# design(dds) <- ~ group
#     #run deseq on the DESeqDataSet results tables for the hippocampus samples
# dds <- DESeq(dds)
# resultsNames(dds)
HIP_5CWT <- DESeq2_restab(c("group", "hippocampus5cKO", "hippocampusWT")    , 3)


############## Get the DEGs #############
source("code/utilities/DESeq2_DEGs.R")
AMY_5CWT_DEGs <- DEGtable2(AMY_5CWT, 4,snakemake@params[["alpha"]])
HIP_5CWT_DEGs <- DEGtable2(HIP_5CWT, 5,snakemake@params[["alpha"]])


############## Gene ontology #############
#make the gene ontology map from enrichR
source('code/utilities/Enrichplot_GO.R')
print(head(AMY_5CWT_DEGs))

#gene ontology log2fc cut off
GO_l2FCcuttoff <- 0.5

ggsave(snakemake@output[[6]], plot = GOmap(AMY_5CWT_DEGs, GO_l2FCcuttoff), width = 9, height = 9)
ggsave(snakemake@output[[7]], plot = GOmap(HIP_5CWT_DEGs, GO_l2FCcuttoff), width = 9, height = 9)

ggsave(snakemake@output[[8]], plot = GOdot(AMY_5CWT_DEGs, GO_l2FCcuttoff), width = 5.5, height = 7)
ggsave(snakemake@output[[9]], plot = GOdot(HIP_5CWT_DEGs, GO_l2FCcuttoff), width = 5.5, height = 7)
