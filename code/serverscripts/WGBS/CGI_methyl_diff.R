# 24.08.08 DNA methylation analysis of bismark bam files at CGI using methylKit
#https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
#https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html
#run from command line by 'Rscript [nameofscript].R'

library(methylKit)
library(ggplot2)
library(genomation)


#we want to compare WT and 5CKO for ESCs and EpiLCs
cells <- c("esc", "EpiLC")


#read in the methylation differences results from the bam files using methRead
#   #samples in list format
#   samp_list <- list()
#   for (s in 1:length(sampCell)){
#     samp_list[[s]] <- sampCell[s]
#   }



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
  
  #combine together, detrand because its cpg
  meth <- unite(myobj, destrand = TRUE)
  print(head(meth))
  
  PCASamples(meth)
  p <- PCASamples(meth)
  
  pdf(paste0("../results/methylKit/PCA_WGBS_", i, ".pdf"), p, width = 4, height = 4)
  dev.off()
  
  #calculate methylation differences
  myDiff <- calculateDiffMeth(meth, mc.cores = 8)
  
  #all of the bases
  myDiff10p <- getMethylDiff(myDiff,difference=10,qvalue=0.1)
    #make a bedgraph file of the differences
  bedgraph(myDiff10p, col.name = "meth.diff", file.name = paste0("../results/methylKit/diff_cpg_10p_", i,".bed"))
  
  #just hypo methylated
  myDiff10p.hypo <- getMethylDiff(myDiff,difference=10,qvalue=0.1,type="hypo")
  
  
  #get the singificantly different cpgs that fall within CGIs
  cpg_anot <- readFeatureFlank("../data/raw/CGI_UCSC.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
  print(head(cpg_anot))  
  
  
  
  diffCpGann <- annotateWithFeatureFlank(as(myDiff10p.hypo,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

  plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")
  q <- plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")
  
  pdf(paste0("../results/methylKit/WGBS_", i, "CGI_hypo_summary.pdf"), p, width = 4, height = 4)
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
  myDiff_islands_10p <- getMethylDiff(myDiff_islands,difference=10,qvalue=0.1)
  
  write.table(myDiff_islands_10p, paste0("../results/methylKit/WGBS_CGI_10_", i, ".csv"), sep = ',', row.names = FALSE, quote = FALSE)
  

  
}

