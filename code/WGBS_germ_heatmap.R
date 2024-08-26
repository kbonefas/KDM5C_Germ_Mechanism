# 2024.08.23 - heatmap of germline gene average % methylation

#want a heatmap of germline genes, split by KDM5C binding or color code?
#TSS heatmap split by CGI-containing vs not

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
plotTSS <- plotTSS[!duplicated(plotTSS$SYMBOL),]
print(paste("rows in merged all germ -", nrow(plotTSS)))
print(head(plotTSS))

write.table(plotTSS, snakemake@output[[1]], sep = ",", row.names = FALSE)


### subset for just genes that have a significant increase CpGme in WT ESC to exEpiLCs
WT_ESCvsEpiLC <- read.csv(snakemake@input[[3]], sep = ",")
WT_ESCvsEpiLC <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$meth.diff > 25,]
print("WT cells")
print(head(WT_ESCvsEpiLC))

#do average between the two samples
#split into groups
#plot just the KDM5C-bound genes?
#make the volcano plot 

### complex heatmap


#need a dataframe that has the row.names as the gene symbols (must be unique) and the % methylation for each sample as columns
hclust_matrix <- data.frame(row.names = plotTSS$SYMBOL, plotTSS[,(ncol(plotTSS)-7):ncol(plotTSS)]) 

#get average of replicates
genotypes <- c("KO", "WT")
cells <- c("esc", "96EpiLC")

avg_hclust_matrix <- data.frame(row.names = plotTSS$SYMBOL, "avg_WT_esc" = rowMeans(subset(hclust_matrix, select = c("WT_A_esc", "WT_B_esc"))), 
"avg_KO_esc" = rowMeans(subset(hclust_matrix, select = c("KO_A_esc", "KO_B_esc"))),
"avg_WT_96EpiLC" = rowMeans(subset(hclust_matrix, select = c("WT_A_96EpiLC", "WT_B_96EpiLC"))),
"avg_KO_96EpiLC" = rowMeans(subset(hclust_matrix, select = c("KO_A_96EpiLC", "KO_B_96EpiLC")))    )


# for(n in 1:4){
# 	for (g in genotypes){
# 		col <- colnames(hclust_matrix)
# 		#get the genotype
# 		col2 <- col[grepl(g, col, fixed = TRUE)]
# 	for (i in cells){
# 		col3 <- col2[grepl(i, col2, fixed = TRUE)]

# 		df <- data.frame(rowMeans(hclust_matrix[,col3]))
# 		colnames(df) <- paste0("avg_", g, "_", i)
# 		print(head(df))
# 		avg[[n]] <- df
		
# 	}}}

print(nrow(avg_hclust_matrix))
print(head(avg_hclust_matrix))

#subset for only the sig increased WT genes
avg_hclust_matrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% WT_ESCvsEpiLC$SYMBOL,]
print("avg hclust_matrix")
print(head(avg_hclust_matrix))

#split by CGI presence
library(ComplexHeatmap)
library(circlize)

heat_region <- function(matrix, CGIstatus){
	#col_fun = colorRamp2(c(0, 50, 100), c("blue", "white", "red"))
	genes <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$Promo_CGI == CGIstatus,]
	print(head(genes))
	plotmatrix <- matrix[row.names(matrix) %in% genes$SYMBOL,]
	col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
	p <- Heatmap(plotmatrix, show_row_names = TRUE, split = 5, column_title = paste("Germline genes with", CGIstatus), column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
	draw(p)
}


#set the column order
#only keep ones that significantly gain me in WT ESC to EpiLCs
#split by CGI and no CGI

#show_row_names = FALSE,
#save the heatmap
pdf(file = snakemake@output[[2]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 30) # The height of the plot in inches
heat_region(avg_hclust_matrix,"no")
dev.off()

pdf(file = snakemake@output[[3]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 20) # The height of the plot in inches
heat_region(avg_hclust_matrix,"CGI")
dev.off()

pdf(file = snakemake@output[[4]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 15) # The height of the plot in inches

genes <- WT_ESCvsEpiLC[WT_ESCvsEpiLC$KDM5C_binding == "Bound",]
print(head(genes))
plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% genes$SYMBOL,]
col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
p <- Heatmap(plotmatrix, column_title = "KDM5C-bound germline genes",show_row_names = TRUE, split = 4, column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
draw(p)
	
dev.off()


#all germline DEGs
EpiLC_DEGs <- read.csv(snakemake@input[[4]], sep = ",")[,"SYMBOL"]
amy_DEGs <- read.csv(snakemake@input[[5]], sep = ",")[,"SYMBOL"]
hip_DEGs <- read.csv(snakemake@input[[6]], sep = ",")[,"SYMBOL"]

DEGs <- unique(c(hip_DEGs, EpiLC_DEGs, amy_DEGs))

pdf(file = snakemake@output[[5]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% DEGs,]
col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
p <- Heatmap(plotmatrix, column_title = "EpiLC germline DEGs",show_row_names = TRUE, split = 2, column_order = c("avg_WT_esc","avg_KO_esc","avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
draw(p)
	
dev.off()



pdf(file = snakemake@output[[6]],   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plotmatrix <- avg_hclust_matrix[row.names(avg_hclust_matrix) %in% DEGs,]
plotmatrix <- plotmatrix[,c("avg_WT_96EpiLC", "avg_KO_96EpiLC")]
col_fun = colorRamp2(range(plotmatrix), hcl_palette = "Reds", reverse = TRUE)
p <- Heatmap(plotmatrix, column_title = "EpiLC germline DEGs",show_row_names = TRUE, split = 2, column_order = c("avg_WT_96EpiLC", "avg_KO_96EpiLC"), heatmap_legend_param = list(title = "% CpGme"), col = col_fun)
draw(p)
	
dev.off()