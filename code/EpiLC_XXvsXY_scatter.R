###24.06.11 scatter plot of DEGs in males vs females

#male vs female germline DEGs
scatterdf <- read.csv(snakemake@input[[1]], sep = ",")
print(head(scatterdf))


#read in the samples results table
samples <- c("XY5cKO", "XX5cHET", "XX5cKO")

for (i in 1:length(samples)){
	restab <- read.csv(snakemake@input[[i+1]], sep = ",", row.names=1)
	restab$ENSEMBL <- row.names(restab)

	#subset for germ genes log2fc, renaming the column to the sample
	restab_new <- subset(restab, select = c("ENSEMBL", "log2FoldChange"))
	colnames(restab_new) <- c("ENSEMBL", paste0(samples[i],"_L2FC") )
	# print(paste(tissues[i], "restab only ENSEMBL and L2FC: "))
	# print(head(restab_new))
	scatterdf <- merge(scatterdf, restab_new, by = "ENSEMBL")
	#print(head(scatterdf))
}

#replace NAs with 0 
scatterdf[is.na(scatterdf)] <- 0

print("plotting df:")
print(head(scatterdf))


#make the edge of the graph the maximum value the points can be, anything outside the edge gets a different shape
edge <- 5
scatterdf[,3:5][scatterdf[,3:5] > edge] <- edge

#scatterdf$axis_range <- ifelse(max(scatterdf[,3:5]) > edge, "Outside axis", "Within axis")
print("axis range df")
print(head(scatterdf))
print(tail(scatterdf))


XYKOvsXXHET_scatter <- subset(scatterdf, Sample_DEG == "All" | Sample_DEG ==  "XY 5CKO and XX 5CHET")
# print("XY 5CKO vs XX 5CHET")
# print(head(XYKOvsXXHET_scatter))

library("ggplot2")
source("code/utilities/colorpalettes.R")
pdf(file = snakemake@output[[1]],   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 2.5) # The height of the plot in inches

ggplot(XYKOvsXXHET_scatter, aes(x=XY5cKO_L2FC, y=XX5cHET_L2FC, color = Sample_DEG)) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_point(size=2.5, alpha = 0.6) +
  xlim(0,edge) +
  ylim(0,edge) +
  xlab("XY 5CKO Log2FC") +
  ylab("XX 5CHET Log2FC") +
  #scale_shape_manual(values = c(19, 9)) +
  scale_color_manual(values = c("#4697fa", "#de5410")) +
  labs(color="Sample DEG", shape="Axis Range")

dev.off()



XYKOvsXXKO_scatter <- subset(scatterdf, Sample_DEG == "All" | Sample_DEG ==  "XY 5CKO and XX 5CKO")
# print("XY 5CKO vs XX 5CKO")
# print(head(XYKOvsXXKO_scatter))

pdf(file = snakemake@output[[2]],   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 2.5) # The height of the plot in inches

ggplot(XYKOvsXXKO_scatter, aes(x=XY5cKO_L2FC, y=XX5cKO_L2FC, color =Sample_DEG)) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_point(size=2.5, alpha = 0.6) +
  xlim(0,edge) +
  ylim(0,edge) +
  xlab("XY 5CKO Log2FC") +
  ylab("XX 5CKO Log2FC") +
  #scale_shape_manual(values = c(19, 9)) +
  scale_color_manual(values = c("#4697fa", "#f28650")) +
  labs(color="Sample DEG", shape="Axis Range")


dev.off()
