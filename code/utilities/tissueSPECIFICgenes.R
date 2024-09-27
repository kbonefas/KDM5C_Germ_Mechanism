#### 23.08.03 expression of tissue specific genes

# Tissue specific gene volcano plot in 5c adult hip and amy


#WT mouse tissue specific genes from https://pubmed.ncbi.nlm.nih.gov/28646208/ 
library("readxl")
Tissue <- read_excel("data/raw/Li_2017_tissueSPECIFICgenes.xlsx")
Tissue <- data.frame(Tissue)
#keep just the tissue and gene ID
Tissue <- Tissue[,1:2]
colnames(Tissue) <- c("tissueID", "ENSEMBL")


#replace the tissue ID with full name and add in colorpalette
#get color palette for tissues
source('code/utilities/colorpalettes.R')

tissueID <- data.frame(tissue = c("Adrenal Gland","Brain","Forestomach","Heart","Kidney","Liver","Large Intenstine","Lung","Muscle","Ovary","Small Intestine","Spleen", "Stomach","Testis","Thymus","Uterus","Vesicular Gland"),tissueID = tiss_list <- unique(Tissue$tissueID), keyvals = tissue_palette)
library(dplyr)
Tissues <- Tissue %>%
  left_join(tissueID, by = c("tissueID" = "tissueID"))
#get rid of old tissue ID
Tissues <- Tissues[,2:4]

#get the number of tissue-specific genes
tissuelist <- unique(Tissues$tissue)
counting <- c()

for (t in tissuelist){
	count = nrow(subset(Tissues, tissue == t))
	counting[t] <- count
}
countingdf <- data.frame(tissue = tissuelist, count = counting)
names(counting) <- tissuelist
print(countingdf)

#get number of DEGs that fall into each tissue category
# n = place in the input snakefile of the DEG csv file
source("code/utilities/parameters.R") #get cut off of log2fc
tissue_genes <- function(n){
	#first read in the DEG table 
	DEGs <- read.csv(snakemake@input[[n]], sep =",")
	names(DEGs)[names(DEGs) == 'X'] <- 'ENSEMBL'
	print(head(DEGs))

	DEGs_up <- subset(DEGs, log2FoldChange > l2fcco)
	#merge together, keeping all columns. Genes that aren't tissue specific will be NAs
	tissue_DEGs <- merge(Tissues, DEGs_up, by = "ENSEMBL", all=T)
	tissue_DEGs <- subset(tissue_DEGs, ENSEMBL %in% DEGs$ENSEMBL)
	# tissue_DEGs <- merge(tissue_DEGs, DEGs, by = c("ENSEMBL", "SYMBOL", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "Direction"), all=T)
	# print(paste(n, "tissue DEGs:", nrow(tissue_DEGs)))
	# print(head(tissue_DEGs))
	return(tissue_DEGs)

}


#volcano plot function

#n = position in snakefile for input DEG table and output volcano plot
#title = name of the dataset for the graph
#td = tissue-specific DEG 
tissue_volcano <- function(n, title, td, offset){
	#read in the results table
	restab <- read.csv(snakemake@input[[n+offset]], sep =",",  row.names = 1)
	#print(head(restab))
	restab$ENSEMBL <- rownames(restab)
	tissue_DEGs <- subset(td, select = c("ENSEMBL", "SYMBOL", "tissue", "keyvals"))
	restab <- merge(tissue_DEGs, restab, by = "ENSEMBL", all=T)
	# print("res tab post merge:")
	# print(head(restab))

	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$keyvals <- ifelse(restab$log2FoldChange >= l2fcco & restab$padj <= 0.1 & is.na(restab$keyvals), 'darkgray',
                         ifelse(restab$log2FoldChange <= -l2fcco & restab$padj <= 0.1  & is.na(restab$keyvals), 'darkgray',
                                restab$keyvals))
	restab$keyvals[is.na(restab$keyvals)] <- 'gray38'
	#print("restab keyvals gray38:")
	#print(head(restab))

	restab$names <- ifelse(restab$keyvals == 'darkgray', 'DEG', 
                         ifelse(restab$keyvals == 'gray38', 'Non-Significant',
                                restab$tissue))



	#for the figure legend
	keyvals <- restab$keyvals
	names(keyvals) <- restab$names 

	# #for the figure legend
	# names(keyvals) <- ifelse(keyvals %in% tissue_DEGs$keyvals, tissue_DEGs$tissue, 
	# 	ifelse(keyvals == 'gray', 'DEG', 'Non-Significant'))

	#Make the axis range narrower so that there is space to see the data
	neglog10 <- function(x){
	  a = log10(x)
	  return(-a)
	}

	l2fcmaxie <- 3
	keyvals.shape <- ifelse(
  		restab$log2FoldChange>l2fcmaxie | neglog10(restab$padj)>10, 9, 16)
	#if the log2fc or the padj is greater than the maximum axis cut off, set it to the cut off and make the point an outlier shape
	restab$log2FoldChange[restab$log2FoldChange>l2fcmaxie] <- l2fcmaxie
	restab$padj[restab$padj< 1e-10] <- 1e-10

	keyvals.shape[is.na(keyvals.shape)] <- 16
	names(keyvals.shape)[keyvals.shape == 9] <- 'Outlier'
	names(keyvals.shape)[keyvals.shape == 16] <- 'Regular'

	#select label
	labels <- c("D1Pas1", "Zar1", "Cyct")

	library(EnhancedVolcano)
	p <- EnhancedVolcano(restab, lab = restab$SYMBOL,
    	                 x = 'log2FoldChange', y = 'padj',
            	         selectLab = labels,
                	     xlim=c(-2,3), ylim=c(0,10),
    	                 xlab = bquote(~Log[2]~ 'fold change'),
        	             ylab = bquote(~-Log[10]~adjusted~italic(P)),
            	         pCutoff = 0.1, FCcutoff = l2fcco,
    	                 title = paste0("Kdm5c-KO ", title, " mRNA-seq"),
        	             subtitle = " ",
            	         labSize = 4.0,
                	     colAlpha = 4/5,
                    	 colCustom = keyvals,
						 shapeCustom = keyvals.shape,
            	         legendPosition = 'right',
                	     gridlines.major = FALSE,
                    	 gridlines.minor = FALSE,
						 drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = 'black',
                	     legendLabSize = 10, legendIconSize = 3.0)

	ggsave(snakemake@output[[n]], plot = p, width = 7, height = 6)


}


#function for histogram of tissue specific DEGs

#td = tissue-specific DEG list
#n = position in snakefile for output histogram
#bary - y-axis of bar graph
hist_tiss <- function(n, title, td, offset, bary){
	library(ggpubr)
	#get the number of each gene
	tiss_count <- c()
	cnt <- 1
	for(t in unique(Tissues$tissue)){
		numb <- nrow(subset(td, td$tissue == t))
		tiss_count[cnt] <- numb
		cnt <- cnt + 1
	}
	tissuedf <- data.frame(Tissue = unique(Tissues$tissue), DEG_number = tiss_count, DEG_percent = tiss_count/nrow(td))
	print("# Tissue Genes:")
	print(tissuedf)

	p <- ggbarplot(tissuedf, "Tissue", "DEG_number", fill = "Tissue", color = "Tissue", palette = tissue_palette, 
				label = TRUE, lab.pos = "out", lab.col = "black") + rotate_x_text(45)
	pl <- ggpar(p, legend = "none", ylab = "# of DEGs", xlab = FALSE,  ylim = c(0, bary), title = paste0("Kdm5c-KO ", title, " Tissue-Specific DEGs") )
	ggsave(snakemake@output[[n+offset]], plot = pl, width = 4, height = 5)

}


####### MA plot of tissue specific genes in 5cKO brain ########
#n = position of results table in snakefile 
#title = name of the dataset for the graph
#td = tissue-specific DEG list
tissue_MA <- function(n, title, td, offset){
	#read in the results table
	restab <- read.csv(snakemake@input[[n+offset]], sep =",",  row.names = 1)
	restab$ENSEMBL <- rownames(restab)
	
	#Get the tissue specific DEGs
	tissue_DEGs <- subset(td, select = c("ENSEMBL", "SYMBOL", "tissue"))
	#annotates which are tissue specific genes, anything that isn't there will be no label
	restab <- merge(tissue_DEGs, restab, by = "ENSEMBL", all=T)


	#if the gene is tissue specific, the color key is the color for the tissue. Otherwise the gene will be gray
	restab$tissue <- ifelse(restab$log2FoldChange >= l2fcco & restab$padj <= 0.1 & is.na(restab$tissue), "Non-specific",
                         ifelse(restab$log2FoldChange <= -l2fcco & restab$padj <= 0.1  & is.na(restab$tissue), "Non-specific",
                                restab$tissue))

	#any remaining NAs are genes that aren't significant - make these a lighter gray
	restab$tissue[is.na(restab$tissue)] <- "Non-significant" 
	names(restab)[names(restab) == 'tissue'] <- 'Tissue' 
	
	#order the tissues for plotting
	#restab$tissue <- factor(restab$tissue, levels = tissue_palette2)

	#calculate the log2basemean
	restab$log2baseMean <- log2(restab$baseMean)

	#set axis limits so spread of data more visible 
	#if the log2fc is greater than the maximum axis cut off, set it to the cut off and make the point an outlier shape
	l2fcmaxie <- 4 #maximum log2fc value
	restab$Shape <- ifelse(restab$log2FoldChange > l2fcmaxie, "Outside axis range", "Within axis range")

	restab$Shape <- factor(restab$Shape, levels = c("Within axis range", "Outside axis range"))
	restab$log2FoldChange[restab$log2FoldChange > l2fcmaxie] <- l2fcmaxie



	print(paste0("plotting " , nrow(restab), " points"))
	#print(head(restab))


	##generate the MA plot
	library("ggpubr")
	q <- ggscatter(restab, x = "log2baseMean", y = "log2FoldChange", color = "Tissue", shape = "Shape", palette = tissue_palette2, title = paste0("Kdm5c-KO ", title, " mRNA-seq"))
		#label = "SYMBOL", repel = TRUE, label.rectangle = TRUE, label.select = c("D1Pas1", "Zar1", "Cyct"))
	q <- q + geom_hline(yintercept=l2fcco, linetype="dashed", color = "gray")
	q <- q + geom_hline(yintercept=-l2fcco, linetype="dashed", color = "gray")
	ggsave(snakemake@output[[n+offset*2]], plot = ggpar(q, legend = "right", ylim = c(-2,l2fcmaxie)), width = 5.5, height = 4)
}



### put all the plotting functions together
#n = position in snakefile for input DEG table and output volcano plot
#title = name of the dataset for the graph
#td = tissue-specific DEG 
#offset - number of samples to shift the input and output numbers to align correctly
tissue_plots <- function(n, title, td, offset, bary){
	#plot the volcano first
	tissue_volcano(n, title, td, offset)
	hist_tiss(n, title, td, offset, bary)
	tissue_MA(n, title, td, offset)

}

### calculate if the tissue is enriched
# td - tissue DEG list
# all tissue genes
# n - position of DEGs
tissue_fisher <- function(td, n, title){
	#get list of DEGs
	DEGs <- read.csv(snakemake@input[[n]], sep =",")
	DEGs <- subset(DEGs, log2FoldChange > l2fcco)
	names(DEGs)[names(DEGs) == 'X'] <- 'ENSEMBL'

	#for every tissue, count how many are DEGs
	pvalues <- c()
	DEGnumb <- c()
	oddsratio <- c()
	cnt <- 1
	for (t in unique(Tissues$tissue)){
		Tissues_genes <- subset(Tissues, tissue == t)
		
		tissue_yes_DEG_yes <- nrow(subset(DEGs, ENSEMBL %in% Tissues_genes$ENSEMBL))

		tissue_no_DEG_no <-  52452 - nrow(Tissues_genes) -  nrow(DEGs) # 52452 is the total number of genes tested

		tissue_no_DEG_yes <- nrow(subset(DEGs, !(ENSEMBL %in% Tissues_genes$ENSEMBL)))

		tissue_yes_DEG_no <- nrow(subset(Tissues_genes, !(ENSEMBL %in% DEGs$ENSEMBL)))

		fisherdf <- data.frame("tissue_yes" = c(tissue_yes_DEG_yes, tissue_yes_DEG_no), "tissue_no" = c(tissue_no_DEG_yes, tissue_no_DEG_no), row.names = c("DEG_yes", "DEG_no"))
		
		fishtest <- fisher.test(fisherdf)
		# print(paste0(title, " fisher results for ", t))
		# print(fisherdf)
		# print(fishtest)

		pvalues[cnt] <- fishtest$p.value
		DEGnumb[cnt] <- tissue_yes_DEG_yes
		oddsratio[cnt] <- fishtest$estimate[[1]]

		cnt <- cnt + 1

	}
		print(paste0(title, " fisher results"))
		fishp <- data.frame(tissue = unique(Tissues$tissue), pvalues = pvalues, DEGnumb = DEGnumb, oddsratio = oddsratio)
		names(fishp)[names(fishp) == 'pvalues'] <- paste0(title, '_pvalues')
		names(fishp)[names(fishp) == 'DEGnumb'] <- paste0(title, '_DEGnumb')
		names(fishp)[names(fishp) == 'oddsratio'] <- paste0(title, '_oddsratio')
		print(fishp)
		return(fishp)

}


