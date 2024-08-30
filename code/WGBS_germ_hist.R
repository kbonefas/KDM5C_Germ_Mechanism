# 2024.08.30 - histogram of germline gene average % methylation

#want a histogram of germline gene methylation, split by CGI vs not, WT vs KO

#read in the KDM5C bound genes/CGI
germ <- read.csv(snakemake@input[[1]], sep = ",")

print(head(germ))

#get the gene coordinate information
source("code/utilities/GeneTSSandTES.R")
germTSS <- geneTSS_df(germ$ENSEMBL, 500)
print(head(germTSS))


## TSS
TSSdf <- read.csv(snakemake@input[[2]], sep = ",")

#rename chr to seqnames to match granges for merging later
names(TSSdf)[names(TSSdf) == "chr"]<- "seqnames"
print(paste("rows in all germ -", nrow(TSSdf)))
print(head(TSSdf))
	
#get the germline genes for each coordinate
plotTSS <- merge(germTSS, TSSdf)
plotTSS <- merge(germ, plotTSS, by = "ENSEMBL")
#remove any duplicated gene names
# plotTSS <- plotTSS[!duplicated(plotTSS$SYMBOL),]
print(paste("rows in merged all germ -", nrow(plotTSS)))
print(head(plotTSS))

write.table(plotTSS, snakemake@output[[1]], sep = ",", row.names = FALSE)


#a dataframe that has the gene symbols and the % methylation for each sample as columns
percmeth <- data.frame(ENSEMBL = plotTSS$ENSEMBL, plotTSS[,(ncol(plotTSS)-7):ncol(plotTSS)]) 

avg_percmeth <- data.frame(ENSEMBL = plotTSS$ENSEMBL, "avg_WT_esc" = rowMeans(subset(percmeth, select = c("WT_A_esc", "WT_B_esc"))), 
"avg_KO_esc" = rowMeans(subset(percmeth, select = c("KO_A_esc", "KO_B_esc"))),
"avg_WT_96EpiLC" = rowMeans(subset(percmeth, select = c("WT_A_96EpiLC", "WT_B_96EpiLC"))),
"avg_KO_96EpiLC" = rowMeans(subset(percmeth, select = c("KO_A_96EpiLC", "KO_B_96EpiLC"))))


### subset for just genes that have a significant increase CpGme in WT ESC to exEpiLCs
WT_ESCvsEpiLC <- read.csv(snakemake@input[[3]], sep = ",")
WT_ESCvsEpiLC <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$meth.diff > 25,]
print(paste0("WT cells: ", nrow(WT_ESCvsEpiLC)))
print(head(WT_ESCvsEpiLC))

#subset for only the sig increased WT genes
avg_percmeth <- avg_percmeth[avg_percmeth$ENSEMBL %in% WT_ESCvsEpiLC$ENSEMBL,]
print("avg_percmeth")
print(head(avg_percmeth))


### want to compare WT vs 5CKO for CGI genes and non-CGI genes
CGI_status <- c("no", "CGI")
library('ggpubr')
source('code/utilities/colorpalettes.R')

plotlist <- list()
stats <- list()


for(s in 1:length(CGI_status)){
	genes <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$Promo_CGI == CGI_status[[s]],]

	avg_percmeth_p <- avg_percmeth[avg_percmeth$ENSEMBL %in% genes$ENSEMBL,]

	plot_df <- data.frame(ENSEMBL = rep(avg_percmeth_p$ENSEMBL, 2), genotype = c(rep("WT", nrow(avg_percmeth_p)), rep("KO", nrow(avg_percmeth_p))), AVGpercentMeth = c(avg_percmeth_p$avg_WT_96EpiLC, avg_percmeth_p$avg_KO_96EpiLC))

	plotlist[[s]] <- gghistogram(plot_df, x = "AVGpercentMeth", add = "mean", rug = TRUE, color = "genotype", fill = "genotype", add_density = TRUE, title = paste("TSS -", CGI_status[s]), palette = c(EpiLC96_XY_KO, EpiLC96_XY_WT))


	#do wt vs KO statistics - Kolmogorov-Smirnov test
	stats[[s]] <- ks.test(plot_df[plot_df$genotype == "WT",][['AVGpercentMeth']], plot_df[plot_df$genotype == "KO",][['AVGpercentMeth']])




}


library("gridExtra")

ggsave(snakemake@output[[2]], grid.arrange(grobs = plotlist, nrow = 1), width = 7, height = 3)

print("no CGI")
print(stats[[1]])

print("CGI")
print(stats[[2]])


# # need to change the dataframe to one column with the sample, one column with the percent methylation, one column with the gene?
# plot_avg_CGI <- data.frame(ENSEMBL = rep(avg_percmeth_CGI$ENSEMBL, 2), genotype = c(rep("WT", nrow(avg_percmeth_CGI)), rep("KO", nrow(avg_percmeth_CGI))), AVGpercentMeth = c(avg_percmeth_CGI$avg_WT_96EpiLC, avg_percmeth_CGI$avg_KO_96EpiLC))
# print("plot_avg_CGI")
# print(head(plot_avg_CGI))
# print(tail(plot_avg_CGI))

# #histogram
# library('ggpubr')
# source('code/utilities/colorpalettes.R')
# p <- gghistogram(plot_avg_CGI, x = "AVGpercentMeth", add = "mean", rug = TRUE, color = "genotype", fill = "genotype", add_density = TRUE, palette = c(EpiLC96_XY_KO, EpiLC96_XY_WT))

# ggsave(snakemake@output[[2]], p, width = 5, height = 3)



# #split by CGI presence
# library(ComplexHeatmap)
# library(circlize)

# heat_region <- function(matrix, CGIstatus){
# 	#col_fun = colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
# 	genes <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$Promo_CGI == CGIstatus,]
# 	print(head(genes))
# 	plotmatrix <- matrix[row.names(matrix) %in% genes$SYMBOL,]
# 	col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
# 	p <- Heatmap(plotmatrix, show_row_names = TRUE, split = 5, column_title = paste("Germline genes with", CGIstatus), column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
# 	draw(p)
# }


# #set the column order
# #only keep ones that significantly gain me in WT ESC to EpiLCs
# #split by CGI and no CGI

# #show_row_names = FALSE,
# #save the heatmap
# pdf(file = snakemake@output[[2]],   # The directory you want to save the file in
#     width = 5, # The width of the plot in inches
#     height = 30) # The height of the plot in inches
# heat_region(avg_hclust_matrix,"no")
# dev.off()

# pdf(file = snakemake@output[[3]],   # The directory you want to save the file in
#     width = 5, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# heat_region(avg_hclust_matrix,"CGI")
# dev.off()

# pdf(file = snakemake@output[[4]],   # The directory you want to save the file in
#     width = 5, # The width of the plot in inches
#     height = 15) # The height of the plot in inches

# genes <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$KDM5C_binding == "Bound",]
# print(head(genes))
# plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% genes$SYMBOL,]
# col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
# p <- Heatmap(plotmatrix, column_title = "KDM5C-bound germline genes",show_row_names = TRUE, split = 4, column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
# draw(p)
	
# dev.off()


# #all germline DEGs
# EpiLC_DEGs <- read.csv(snakemake@input[[4]], sep = ",")[,"SYMBOL"]
# amy_DEGs <- read.csv(snakemake@input[[5]], sep = ",")[,"SYMBOL"]
# hip_DEGs <- read.csv(snakemake@input[[6]], sep = ",")[,"SYMBOL"]

# DEGs <- unique(c(hip_DEGs, EpiLC_DEGs, amy_DEGs))

# pdf(file = snakemake@output[[5]],   # The directory you want to save the file in
#     width = 5, # The width of the plot in inches
#     height = 10) # The height of the plot in inches

# plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% DEGs,]
# col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
# p <- Heatmap(plotmatrix, column_title = "EpiLC germline DEGs",show_row_names = TRUE, split = 2, column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
# draw(p)
	
# dev.off()



# pdf(file = snakemake@output[[6]],   # The directory you want to save the file in
#     width = 5, # The width of the plot in inches
#     height = 10) # The height of the plot in inches

# plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% DEGs,]
# plotmatrix <- plotmatrix[,c("avg_WT_96EpiLC", "avg_KO_96EpiLC")]
# col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
# p <- Heatmap(plotmatrix, column_title = "EpiLC germline DEGs",show_row_names = TRUE, split = 2, column_order = c("avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
# draw(p)
	
# dev.off()