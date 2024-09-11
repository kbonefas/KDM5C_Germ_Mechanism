# 24.08.24 DNA methylation analysis of bismark bam files using methylKit
    #https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
#run from command line by 'Rscript [nameofscript].R'

library(methylKit)
library(ggplot2)
# library(stringr)
library(ChIPseeker)
library(dplyr)

#need granges as input for the bedfile
#read in the bed file of CGIs
bed <- readPeakFile("../data/raw/CGI_UCSC.bed")
print(head(bed))


#read in the bam files
bamfiles <- list.files(path = "../data/processed/BAM_dedup", pattern = "*.sorted.dedup.bam",
           full.names = TRUE)

#remove the bam.bai files
bamfiles <- bamfiles[!endsWith(bamfiles, '.bai')]


#we want to compare WT and 5CKO for ESCs and EpiLCs
cells <- c("esc", "EpiLC")
#empty lists to store the output

mylists <- list(list(), list())
names(mylists) <- cells

for (i in cells){
    #subset for samples of that cell type
    sampCell <- grep(i, bamfiles, value=TRUE)

    #for every sample in sampcell, get the treatment (WT or KO)
    #get the sample ID
    ID <- gsub("../data/processed/BAM_dedup/", "", sampCell)
    ID <- gsub(".sorted.dedup.bam", "", ID)
    
    #format as a list
    ID_list <- list()
    for (d in 1:length(ID)){
        ID_list[[d]] <- ID[d]
    }
    # print("IDs as list")
    # print(ID_list)
    
    #get the sample genotype (first two characters) - vector format
    geno <- substr(ID, start = 1, stop = 2)
    # print("treatments")
    # print(treat)
    # treatment is 1 or 0, 1 for KO, 0 for WT
    treat <- c()
    for (g in 1:length(geno)) {
        treat[g] <- ifelse(geno[g] == "KO", 1, ifelse(geno[g] == "WT", 0, NA))
    }
    print(treat)

    #samples in list format
    samp_list <- list()
    for (s in 1:length(sampCell)){
        samp_list[[s]] <- sampCell[s]
    }


    #save the alignment into the correct list
    #When it's run and if you save the results you can then read the results in directly using methRead
    mylists[[i]] <- processBismarkAln(location = samp_list, nolap = TRUE, mincov = 3,
		                    sample.id = ID_list, assembly="mm10", treatment = treat, 
                          read.context="CpG", save.folder="../results/methylKit")
    

}



for (n in cells){
	# #get the methylation for cgis
	# regions <- regionCounts(mylists[[n]],bed)

	# united_regions <- unite(regions, destrand=FALSE)
	# print("united regions")
	# print(head(united_regions))

	# write.table(united_regions, paste0("../results/methylKit/regionCounts_CGI_3min_", n, ".csv"), sep = ',', row.names = FALSE, quote = FALSE)


    #combine methylation calls into one
    meth <- unite(mylists[[n]], destrand = TRUE)
    print(head(meth))

    #get the PCA plot
    PCASamples(meth)
	 
	pooled.meth <- pool(meth, sample.ids=c("KO","WT"))
	dm.pooledf <- calculateDiffMeth(pooled.meth, mc.cores = 8)

	myDiff10p <- getMethylDiff(dm.pooledf, difference=10, qvalue=0.01)
    #make a bedgraph file of the differences
	bedgraph(myDiff10p, col.name = "meth.diff", file.name = paste0("../results/methylKit/bedgraph_diff_cpg_10p_q01_3min_pooled_", n,".bed"))
  
    write.table(myDiff10p, paste0("../results/methylKit/WGBS_CpG_10p_q01_3min_", n, ".csv"), sep = ',', row.names = FALSE, quote = FALSE)

}

