#23.11.02 Analysis of germline DEG expression in testis cells

############## mueller wwv testis ###############
mueller <- read.csv(snakemake@input[[1]], header=TRUE,sep="\t", row.names=1)
print(head(mueller))

#get the average of WT and WWv into a new data frame
fpkm_WT <- apply(mueller[1:2], 1, mean)
fpkm_WWv <- apply(mueller[3:4], 1, mean)
avg.fpkm <- cbind(fpkm_WT, fpkm_WWv)
colnames(avg.fpkm) <- c("WT","WWv")
avg.fpkm <- as.data.frame(avg.fpkm)

#get the germline DEGs
samplelist <- c("amy", "hip")
testisDEGs <- read.csv(snakemake@input[[3]]) 

print(head(testisDEGs))
print(nrow(testisDEGs))

#get their expression in WT vs WWv
wwv_fpkm <- subset(avg.fpkm, rownames(avg.fpkm) %in% testisDEGs$ENSEMBL)

#to make the points more visible, set outlier values to a different shape
wwv_fpkm$shape <- ifelse(wwv_fpkm$WT>300, 9, 16)
#set any values greater than 300 to 300
wwv_fpkm$WT[wwv_fpkm$WT>300] <- 300

library("ggplot2")

pdf(file = snakemake@output[[1]],   # The directory you want to save the file in
    width = 4.25, # The width of the plot in inches
    height = 4) # The height of the plot in inches
ggplot(wwv_fpkm, aes(x=WT, y=WWv)) + 
  geom_point(size = 3, color = "#2A7AFF", shape = wwv_fpkm$shape, show.legend = TRUE)+
  geom_abline(intercept = 0, slope = 1, 
              linetype="dashed", size=0.75)+
  xlim(c(0, 300)) + ylim(c(0,300))+
  labs(x = "Wild-type FPKM", y = "W/Wv (Germ cell depleted) FPKM", title = "Expression of testis DEGs \n with germ cell depletion" )
dev.off()



############# green 2018 scRNAseq testis ##############
green <- read.csv(snakemake@input[[2]], header=TRUE,sep="\t", row.names=NULL)
print(head(green))

#replace some names for cleaner formatting
green[green == "Elongating"] <- "Elongating Spermatid"
green[green == "RoundSpermatid"] <- "Round Spermatid"
green[green == "Unknown"] <- "Mesenchymal" #identified in the paper
green[green == "InnateLymphoid"] <- "Innate Lymphoid"


celltype_count <- c()

for (c in unique(green$CellType)){
	greengenes <- subset(green, CellType == c)
	count <- nrow(subset(testisDEGs, SYMBOL %in% greengenes$gene))
	celltype_count <- append(celltype_count, count) 
}

celltype_df <- data.frame(cell_type = unique(green$CellType), count = celltype_count)
print(celltype_df)

library('ggpubr')
p <- ggbarplot(celltype_df, x = "cell_type", y = "count", color = "cell_type", fill = "cell_type", label = TRUE, lab.pos = "out", title = "Expression of testis DEGs in testis cell types")
ggsave(snakemake@output[[2]], ggpar(p, x.text.angle = 45, legend = "none", xlab = FALSE, ylab = "# of testis DEGs" ), width = 5, height = 5)