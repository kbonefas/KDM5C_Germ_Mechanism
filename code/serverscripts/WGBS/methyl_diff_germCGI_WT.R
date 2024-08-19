# 24.08.08 DNA methylation analysis of bismark bam files at CGI using methylKit
#https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
#https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html
#https://compgenomr.github.io/book/extracting-interesting-regions-differential-methylation-and-segmentation.html
#run from command line by 'Rscript [nameofscript].R'

#qvalue cut off
qval <- 0.01


library(methylKit)
library(ggplot2)
library(genomation)



#regional summary analysis
library(ChIPseeker)
library(dplyr)

#need granges as input for the bedfile
#read in the bed file of CGIs for germline gene promoters
bed <- readPeakFile("../data/raw/CGI_UCSC_all_germ.bed")
print(head(bed))


#get list of WT files
difffiles <- Sys.glob("../results/methylKit/WT*_CpG.txt")
difffiles
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


#get the sample cell type (first two characters) - vector format
cell <- rep(c("EpiLC", "ESC"), 2)

# treatment is 1 or 0, 1 for KO, 0 for WT
treat <- c()
for (g in 1:length(cell)) {
treat[g] <- ifelse(cell[g] == "EpiLC", 1, ifelse(cell[g] == "ESC", 0, NA))
}
  
  #read in the files
myobj <- methRead(location = samp_list, sample.id = ID_list, assembly="mm10", treatment = treat,
                   context="CpG", mincov = 3)
#filter based on read coverage
myobj <- filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                    hi.count=NULL,hi.perc=99.9)
  

#get the methylation for cgis
regions <- regionCounts(myobj, bed)
# print("all regions")
# print(head(regions))

united_regions <- unite(regions, destrand=FALSE)
  #can try setting min per group to 1 so bases only have to be covered in one sample/group to be counted
print("united regions")
print(head(united_regions))

#then do differential methylation for the regions
dm.regions <- calculateDiffMeth(united_regions, mc.cores = 8)
#this should be like a DESEQ2 results table
head(dm.regions)

#just saving the one methylation results. Want the average methylation
write.table(dm.regions, "../results/methylKit/WGBS_restab_regionCounts_germCGI_WT_min3_ESCvsEpiLC.csv", sep = ',', row.names = FALSE, quote = FALSE)

#get just the significant ones - should be like the DESeq2 DEGs
myDiff_regions <- getMethylDiff(dm.regions, difference=10, qvalue=qval)
write.table(myDiff_regions, "../results/methylKit/WGBS_getmethyldiff_germCGI_p10_min3_q01_WT_ESCvsEpiLC.csv", sep = ',', row.names = FALSE, quote = FALSE)



#pool the samples for calculating the methylation differences
#unite the samples
# meth <- methylKit::unite(myobj, destrand=TRUE)
# pooled.meth <- pool(meth, sample.ids=c("EpiLC","ESC"))
# dm.pooledf <- calculateDiffMeth(pooled.meth, mc.cores = 8)

#methylation differences pooled
united_dmC <- unite(myobj, destrand=FALSE)

dm.C <- calculateDiffMeth(united_dmC, mc.cores = 8)

#all of the bases
myDiff10p <- getMethylDiff(dm.C, difference=10, qvalue=qval)
#make a bedgraph file of the differences
bedgraph(myDiff10p, col.name = "WT ESC vs EpiLC min3 q01 unpooled", file.name = "../results/methylKit/bedgraph_diff_WT_ESCvsEpiLC_p10_min3_q01_notpooled.bed")  
  
# #get the singificantly different cpgs that fall within CGIs
# cpg_anot <- readFeatureFlank("../data/raw/CGI_UCSC.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
# print(head(cpg_anot))  

# myDiff10p_hyper <- getMethylDiff(dm.pooledf, difference=10, qvalue=qval, type = "hyper")

# diffCpGann <- annotateWithFeatureFlank(as(myDiff10p_hyper,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

# plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")
# q <- plotTargetAnnotation(diffCpGann, main = "Differential Methylation Annotation")

# pdf(paste0("../results/methylKit/WGBS_summary_CGI_hyper_ESCvsEpiLC.pdf"), q, width = 4, height = 4)
# dev.off()
  

  
# #get the list of differentially methylated CGIs
# myobj_islands <- regionCounts(myobj, cpg_anot$CpGi)
# # Filter the summarized counts by coverage
# myobj_islands_filt <- filterByCoverage(myobj_islands, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
# # Perform simple normalization 
# myobj_islands_filt_norm <- normalizeCoverage(myobj_islands_filt, method = "median")
# # Merge the samples again
# meth_islands <- unite(myobj_islands_filt_norm, destrand=TRUE)
  
# #differentially methylation for islands
# myDiff_islands <- calculateDiffMeth(meth_islands)
#   # Rank by significance
# myDiff_islands <- myDiff_islands[order(myDiff_islands$qvalue),]
#   # get all differentially methylated CpG Islands
# myDiff_islands_10p <- getMethylDiff(myDiff_islands, difference=10, qvalue=qval)
  
# write.table(myDiff_islands_10p, paste0("../results/methylKit/WGBS_getmethyldiff_germCGI_10_WT_ESCvsEpiLC.csv"), sep = ',', row.names = FALSE, quote = FALSE)
  


