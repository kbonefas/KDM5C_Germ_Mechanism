#############################################################################################
#                                                                                           #
# 2022 12 14 DESeq2 using snakemake of Kdm5c mutant EpiLCs                                  #
# Inputs:
#       cts = counts file (.txt)                                                            #
#       sampleinfo = sample information sheet #MAKE SURE THE SAMPLES MATCH THE CTS FILE NAMES
# Libraries: DESeq2, ggplot2                                                                #
# Outputs: PCA plot of samples and DESeq2 results tables (.csv)                             #
#                                                                                           #
#############################################################################################
#load the color palette utilities first
source('code/utilities/colorpalettes.R')

###1) format the raw counts table so that all values are integers and there are no gene variant numbers
cts <- read.csv(snakemake@input[["cts"]], sep ="\t", row.names = 1)
    #convert counts to integers
cts[1:ncol(cts)] <- lapply(cts[1:ncol(cts)], as.integer)

###2) Make the sample sheet (a csv file that has the genotype and sex for each sample)
#get the column (sample) names
samples <-  colnames(cts)

#get the sex of each sample
sex <- c()
ct <- 1
for(r in samples){
		#get the first two characters of the string. If they're XX, they're female
	if (substr(r, 1, 2) == "XX"){
		rei <- "XX"
	    } else if (substr(r, 1, 2) == "XY"){ #if they're male
		rei <- "XY"
	    } else {
		rei <- "ERROR"
	}
	sex[ct] <- rei
	ct <- ct + 1
	
}

#get the genotype of each sample
genotype <- c()
cnt <- 1
for(g in samples){
		#get the end fo the string
	if (substr(g, 10, 11) == "WT"){
		geno <- "WT"
        } else if (substr(g, 10, 12) == "HET") {
        geno <- "5cHET"
        } else if (substr(g, 10, 11) == "KO"){
        geno <- "5cKO"
        } else { #the rest of the samples have 7 characters
		geno <- "ERROR"
	}
	genotype[cnt] <- geno
	cnt <- cnt + 1
}


SampleInfo <- data.frame(Sample = colnames(cts), Genotype = genotype, Sex = sex, Type = rep("paired.end", length(genotype)))
write.table(SampleInfo, snakemake@output[[1]], sep=",", row.names=FALSE)


###3) Set up and run DESeq2 
library('DESeq2')
    #create the coldata
coldata <- data.frame(row.names = SampleInfo$Sample, genotype = SampleInfo$Genotype, sex = SampleInfo$Sex, type = SampleInfo$Type)
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

    #add in sex so that we can account for the interaction between genotype and sex
dds$group <- factor(paste0(dds$sex, dds$genotype))
    #new design using the group variable
dds$group <- relevel(dds$group, ref = "XYWT")
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
    #adding sex as a separate variable back in so that the shape of the point is dependent upon sex
data <- data.frame(data, Sex = SampleInfo$Sex)
q <- ggplot(data, aes(PC1, PC2, color=group, shape = Sex)) +
    geom_point(size=4.5) +
    ggtitle("PCA Plot of EpiLCs")+
    scale_color_manual(values=EpiLCpalette) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) 

    #Save the PCA plot
ggsave(snakemake@output[[2]], plot = q, width = 5, height = 4)


###5) generate results tables for each comparison  
    #This function runs the results and the shrink for each comparison
DESeq2_restab <- function(n, i){
    res_unshrunk <- results(dds, name=n, alpha = snakemake@params[["alpha"]])
    #LFCshrink to normalize small count variation
    res <- lfcShrink(dds, coef=n, res=res_unshrunk, type = "ashr")
    #prints the results summary (# of DEGs)
    summary(res)
    #i is the number in the output list for the snakefile that corresponds to the samples
    write.table(res, snakemake@output[[i]], sep=",", col.names=NA)
    return(res)
}
    

XY_KOWT <- DESeq2_restab("group_XY5cKO_vs_XYWT", 3)

dds$group <- relevel(dds$group, ref = "XXWT")
design(dds) <- ~ group
    #run deseq on the DESeqDataSet results tables for the female samples
dds <- DESeq(dds)
XX_HETWT <- DESeq2_restab("group_XX5cHET_vs_XXWT", 4)
XX_KOWT <- DESeq2_restab("group_XX5cKO_vs_XXWT", 5)

############# Get the DEGs #############


source("code/utilities/parameters.R") #get cut off of log2fc
source("code/utilities/DESeq2_DEGs.R")
invisible(DEGtable3(XY_KOWT, 6, snakemake@params[["alpha"]], l2fcco_ESCEpiLC))
invisible(DEGtable3(XX_HETWT, 7, snakemake@params[["alpha"]], l2fcco_ESCEpiLC))
invisible(DEGtable3(XX_KOWT, 8, snakemake@params[["alpha"]], l2fcco_ESCEpiLC))
