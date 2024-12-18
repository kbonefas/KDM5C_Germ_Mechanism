# 24.08.15 DNA methylation analysis of bismark bam files at CGI using methylKit
#https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
#https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html
#https://compgenomr.github.io/book/extracting-interesting-regions-differential-methylation-and-segmentation.html
#run from command line by 'Rscript [nameofscript].R'

## calculating methylation differences in WT vs Ko ESCs and EpiLCs

library(methylKit)
library(ggplot2)
library(genomation)



#regional summary analysis
library(ChIPseeker)
library(dplyr)

#need granges as input for the bedfile
#read in the bed file of CGIs
bed <- readPeakFile("../data/raw/CGI_UCSC.bed")
print(head(bed))


#we want to compare WT and 5CKO for ESCs and EpiLCs
cells <- c("esc", "EpiLC")


for (i in cells){
  #get list of files
  difffiles <- list.files(path = "../results/methylKit", pattern = paste0("*",i,"_CpG.txt"),
                         full.names = TRUE)
  #make as a list
  samp_list <- list()
  for (s in 1:length(difffiles)){
    samp_list[[s]] <- difffiles[s]
  }
  
  #get the sample names (IDs)   
  ID <- gsub("../results/methylKit/", "", difffiles)
  ID <- gsub("_CpG.txt", "", ID)
  #format as a list
  ID_list <- list()
  for (d in 1:length(ID)){
    ID_list[[d]] <- ID[d]
  }

  print("sample ID_list")
  print(ID_list)


  #get the sample genotype (first two characters) - vector format
  geno <- substr(ID, start = 1, stop = 2)

  # treatment is 1 or 0, 1 for KO, 0 for WT
  treat <- c()
  for (g in 1:length(geno)) {
    treat[g] <- ifelse(geno[g] == "KO", 1, ifelse(geno[g] == "WT", 0, NA))
  }
  
  #read in the files
  myobj <- methRead(location = samp_list, sample.id = ID_list, assembly="mm10", treatment = treat,
                   context="CpG", mincov = 5)
  #filter based on read coverage
  myobj <- filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
  

  #get the methylation for cgis
  regions <- regionCounts(myobj,bed)
  print(head(regions))


  united_regions <- unite(regions, destrand=FALSE)
	print("united regions")
	print(head(united_regions))

	write.table(united_regions, paste0("../results/methylKit/regionCounts_CGI_5minCOV_", i, "_WTvsKO.csv"), sep = ',', row.names = FALSE, quote = FALSE)


  #pool the samples for calculating the methylation differences
    #unite the samples
  meth <- methylKit::unite(myobj, destrand=TRUE)
  pooled.meth <- pool(meth, sample.ids=c("KO","WT"))
  dm.pooledf <- calculateDiffMeth(pooled.meth, mc.cores = 8)

  #all of the bases
  myDiff10p <- getMethylDiff(dm.pooledf, difference=10, qvalue=0.1)
    #make a bedgraph file of the differences
  bedgraph(myDiff10p, col.name = "meth.diff", file.name = paste0("../results/methylKit/diff_cpg_10p_q1_pooled_", i,"_WTvsKO.bed"))
  
  
  #get the singificantly different cpgs that fall within CGIs
  cpg_anot <- readFeatureFlank("../data/raw/CGI_UCSC.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
  print(head(cpg_anot))  
    
  
  diffCpGann <- annotateWithFeatureFlank(as(myDiff10p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

  plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")
  q <- plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")
  
  pdf(paste0("../results/methylKit/WGBS_", i, "_WTvsKO_CGI_summary.pdf"), q, width = 4, height = 4)
  dev.off()
  

  #get the list of differentially methylated CGIs
  myobj_islands <- regionCounts(myobj, cpg_anot$CpGi)
  # Filter the summarized counts by coverage
  myobj_islands_filt <- filterByCoverage(myobj_islands, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
  # Perform simple normalization 
  myobj_islands_filt_norm <- normalizeCoverage(myobj_islands_filt, method = "median")
  # Merge the samples again
  meth_islands <- unite(myobj_islands_filt_norm, destrand=TRUE)
  
  #differentially methylation for islands
  myDiff_islands <- calculateDiffMeth(meth_islands)
  # Rank by significance
  myDiff_islands <- myDiff_islands[order(myDiff_islands$qvalue),]
  # get all differentially methylated CpG Islands
  myDiff_islands_10p <- getMethylDiff(myDiff_islands, difference=10, qvalue=0.1)
  
  write.table(myDiff_islands_10p, paste0("../results/methylKit/WGBS_CGI_10p_q1", i, "_WTvsKO.csv"), sep = ',', row.names = FALSE, quote = FALSE)
  
  
}






