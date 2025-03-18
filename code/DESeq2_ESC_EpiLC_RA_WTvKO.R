#2310039 Analysis of ESC to EpiLC RNAseq

#questions we want to address:
	#what germline genes are differentially expressed in 5CKO ESCs and EpiLCs at each time point
	#are those genes still expressed without RA?
		#plot log2scale of counts/TPM? Heatmap of scaled TPM? Example gene?


##1) read in and format the counts table
cts_all <- read.csv(snakemake@input[["cts"]], sep ="\t", row.names = 1)
    #remove unused columns with gene names
cts <- cts_all[4:ncol(cts_all)] 
print("initial cts")
print(head(cts))

##2) read in the sample information sheet (CSV file)
SampleInfo <- read.csv(snakemake@input[["sampleinfo"]], sep =",") 
print(SampleInfo)
##3) Set up and run DESeq2

#create the coldata
coldata <- data.frame(row.names = SampleInfo$ID, genotype = factor(SampleInfo$Genotype), timepoint = factor(SampleInfo$Timepoint), RA = SampleInfo$RA, batch = factor(SampleInfo$Batch))
print(coldata)

#order the column data to match the order of the samples in the cts file
	#first check if sample names are in the column names of the counts
print(all(rownames(coldata) %in% colnames(cts)))
	#if true, subset for cts columns within the rownames of coldata
cts <- cts[, rownames(coldata)]
print(head(cts))
	#check that they match
print(all(rownames(coldata) == colnames(cts)))

library('DESeq2')

#run deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                        colData = coldata,
                        design = ~ batch + genotype)

#add in the other variables
dds$group <- factor(paste0(dds$genotype, dds$timepoint, dds$RA))
design(dds) <- ~ group
#dds$genotype <- relevel(dds$genotype, ref ="WT")
dds <- DESeq(dds)

##4) make a PCA plot 
#variance stabilizing transformation
    #essentially accounts for sampling variability of low counts while also normalizing variability across all samples 
vsd <- vst(dds, blind=TRUE)
data <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
   #adding region as a separate variable back in so that the shape of the point is dependent upon RA
data <- data.frame(data, RA = coldata$RA, genotime = paste0(coldata$genotype, coldata$timepoint))
print(head(data))

data$RA[data$RA == 'RA'] <- 'RA +'
data$RA[data$RA == 'NO'] <- 'RA -'
data$RA[data$RA == ''] <- 'ESC'

library('ggplot2')
library('ggpubr')

source("code/utilities/colorpalettes.R")
q <- ggscatter(data, "PC1", "PC2", color = "genotime", size = 4.5, palette = allRA, shape = "RA", title = "PCA of XY ESCs to EpiLCs, +/- Retinoic Acid", xlab = paste0("PC1: ",percentVar[1],"% variance"), ylab = paste0("PC2: ",percentVar[2],"% variance"))
q <- q + theme_classic()
	#Save the PCA plot
ggsave(snakemake@output[[1]], plot = ggpar(q, legend = "right"), width = 5, height = 4)


##5) Generate the results tables for the contrasts of interest
germ <- read.csv(snakemake@input[["germ"]], sep = ",")
source("code/utilities/parameters.R")
l2fcco <- l2fcco_ESCEpiLC

ratimes <- c("0", "48RA", "48NO", "96RA", "96NO")
germDEGlist <- list()
allgermDEGs <- data.frame()
#all the results tables dataframes
resdfs <- list()
for (i in 1:length(ratimes)){
	#timepoints you're comparing
	KOsample <- paste0("5CKO",ratimes[i])
	WTsample <- paste0("WT", ratimes[i])
	res_unshrunk <- results(dds, contrast=c("group", KOsample, WTsample), alpha = snakemake@params[["alpha"]])
	res <- lfcShrink(dds, contrast=c("group", KOsample, WTsample), res=res_unshrunk, type = "ashr")
	print(paste("for", ratimes[i]))
	print(summary(res))
	write.table(res, snakemake@output[[i+1]], sep=",", col.names=NA)

	#next get the germline DEGs
	res <- data.frame(res)
	resdfs[[i]] <- res
	germDEGs <- subset(res, row.names(res) %in% germ$ENSEMBL & res$padj < snakemake@params[["alpha"]] & log2FoldChange > l2fcco)
	germDEGs$ENSEMBL <- row.names(germDEGs)
	germDEGs <- merge(germ, germDEGs, by = "ENSEMBL")
	print(paste("There are", nrow(germDEGs), "germline DEGs for", ratimes[i]))
	print(head(germDEGs))
	write.table(germDEGs, snakemake@output[[i+6]], sep=",", row.names =FALSE)
	germDEGlist[[i]] <- germDEGs

	#add the ENSEMBL and SYMBOl information to the germline DEG list
	allgermDEGs <- rbind(allgermDEGs, germDEGs[,1:2])
}

ratimesnames <- c("ESC", "RA48", "NO48", "RA96", "NO96")
names(germDEGlist) <- ratimesnames
names(resdfs) <- ratimesnames



