# 24.08.19 DNA methylation analysis of bismark bam files at germline CGI and promoter regions using methylKit
# 5CKO vs WT exEpiLCs
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



#get list of WT files
difffiles <- Sys.glob("../results/methylKit/*EpiLC_CpG.txt")
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
cell <- c("KO", "KO", "WT", "WT")

# treatment is 1 or 0, 1 for KO, 0 for WT
treat <- c()
for (g in 1:length(cell)) {
treat[g] <- ifelse(cell[g] == "KO", 1, ifelse(cell[g] == "WT", 0, NA))
}
  
  #read in the files
myobj <- methRead(location = samp_list, sample.id = ID_list, assembly="mm10", treatment = treat,
                   context="CpG", mincov = 3)
#filter based on read coverage
myobj <- filterByCoverage(myobj,lo.count=3,lo.perc=NULL,
                    hi.count=NULL,hi.perc=99.9)
  

#get the differential methylation for the bedfile regions of interest
#outname is the output file name
METHoverBED <- function(bedfile, outname){
	bed <- readPeakFile(bedfile)
	print(head(bed))

	#count the methylation in the regions for each sample
	regions <- regionCounts(myobj, bed)

	#unite the samples together so only regions covered in all samples are calculated
		#can try setting min per group to 1 so bases only have to be covered in one sample/group to be counted
		#only use destrand with CpGme
	united_regions <- unite(regions, destrand=TRUE)
  	
	#then do differential methylation for the regions
	dm.regions <- calculateDiffMeth(united_regions, mc.cores = 8)
	#this should be like a DESEQ2 results table
	head(dm.regions)

	#methylation results for all regions, similar to results table for DESeq2
	write.table(dm.regions, paste0("../results/methylKit/WGBS_restab_regionCounts_", outname,"_min3_WTvsKO_EpiLC.csv"), sep = ',', row.names = FALSE, quote = FALSE)

	#get just the significant ones - should be like the DESeq2 DEGs
	myDiff_regions <- getMethylDiff(dm.regions, difference=25, qvalue=qval)
	write.table(myDiff_regions, paste0("../results/methylKit/WGBS_getmethyldiff_", outname,"_p25_min3_q01_WTvsKO_EpiLC.csv"), sep = ',', row.names = FALSE, quote = FALSE)

}


#location of bedfiles you want to assess methylation over
bedregions <- c("../data/raw/TSS_window_500bp_all_genes.bed")
names(bedregions) <- c("allgenesTSS500")

for (r in 1:length(bedregions)){
	METHoverBED(bedregions[[r]], names(bedregions)[[r]])
}

# #looking at all the CpG bases to make a bedrgaph file of differentially methylated Cs
# #methylation differences
# united_dmC <- unite(myobj, destrand=TRUE)
# dm.C <- calculateDiffMeth(united_dmC, mc.cores = 8)
# myDiff10p <- getMethylDiff(dm.C, difference=25, qvalue=qval)
# #make a bedgraph file of the differences
# bedgraph(myDiff10p, col.name = "meth.diff", file.name = "../results/methylKit/bedgraph_diff_WTvsKO_EpiLC_p25_min3_q01_notpooled.bed")  
  






#pool the samples for calculating the methylation differences
#unite the samples
# meth <- methylKit::unite(myobj, destrand=TRUE)
# pooled.meth <- pool(meth, sample.ids=c("EpiLC","ESC"))
# dm.pooledf <- calculateDiffMeth(pooled.meth, mc.cores = 8)




