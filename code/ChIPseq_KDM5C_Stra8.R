#24.07.29 - Do kdm5c bound vs unbound correlate with stra8 targest

#Kdm5c bound germ genes
germ_KDM5C <- read.csv(snakemake@input[[1]], sep = ",")

#Stra8 bound genes
Stra8 <- read.csv(snakemake@input[[2]], sep = ",")

#get stra8-bound genes that are either KDM5C bound or unbound
germ_KDM5C_Stra8 <- merge(germ_KDM5C, Stra8, by = "ENSEMBL")


#don't have every germline-enriched gene unfortunately (1130 of 1287)
Stra8df <- data.frame(Kdm5c_binding = c(rep("Bound", 2), rep("Unbound", 2)),  Stra8_binding = rep(c("Bound", "Unbound"), 2))
Stra8df$Stra8_Count <- c(nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Bound" & germ_KDM5C_Stra8$STRA8_bound == "STRA8-bound")) , nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Bound" & germ_KDM5C_Stra8$STRA8_bound == "not STRA8-bound")), nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8$STRA8_bound == "STRA8-bound")), nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8$STRA8_bound == "not STRA8-bound")))

Stra8df$Stra8_Percent_plot <- ifelse(Stra8df$Kdm5c_binding == "Bound", Stra8df$Stra8_Count/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Bound")$Stra8_Count) * 100, ifelse(Stra8df$Kdm5c_binding == "Unbound", Stra8df$Stra8_Count/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Unbound")$Stra8_Count) * 100, 0))




#cpg island near KDM5C-bound promoter
Stra8df$CpG_island <- rep(c("noCGI", "CGI"), 2)
Stra8df$CGI_Count <- c(nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Bound" & germ_KDM5C_Stra8$Promoter_CpG_island == "noCGI")) , nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Bound" & germ_KDM5C_Stra8$Promoter_CpG_island == "CGI")), nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8$Promoter_CpG_island == "noCGI")), nrow(subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8$Promoter_CpG_island == "CGI")))


Stra8df$CGIPercent_plot <- ifelse(Stra8df$Kdm5c_binding == "Bound", Stra8df$CGI_Count/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Bound")$CGI_Count) * 100, ifelse(Stra8df$Kdm5c_binding == "Unbound", Stra8df$CGI_Count/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Unbound")$CGI_Count) * 100, 0))



Stra8df$Stra8_Percent <- as.integer(round(Stra8df$Stra8_Percent_plot))
Stra8df$CGI_Percent <- as.integer(round(Stra8df$CGIPercent_plot))


#do for DEGs
#get the DEG list
amyDEGs <- read.csv(snakemake@input[[3]], sep = ",")$ENSEMBL
hipDEGs <- read.csv(snakemake@input[[4]], sep = ",")$ENSEMBL
EpiLCDEGs <- read.csv(snakemake@input[[5]], sep = ",")$ENSEMBL
allDEGs <- unique(c(amyDEGs, hipDEGs, EpiLCDEGs))

germ_KDM5C_Stra8_DEGs <- subset(germ_KDM5C_Stra8, germ_KDM5C_Stra8$ENSEMBL %in% allDEGs)


Stra8df$Stra8_Count_DEGs <- c(nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Bound" & germ_KDM5C_Stra8_DEGs$STRA8_bound == "STRA8-bound")) , nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Bound" & germ_KDM5C_Stra8_DEGs$STRA8_bound == "not STRA8-bound")), nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8_DEGs$STRA8_bound == "STRA8-bound")), nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8_DEGs$STRA8_bound == "not STRA8-bound")))

Stra8df$Stra8_Percent_DEGs_plot <- ifelse(Stra8df$Kdm5c_binding == "Bound", Stra8df$Stra8_Count_DEGs/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Bound")$Stra8_Count_DEGs) * 100, ifelse(Stra8df$Kdm5c_binding == "Unbound", Stra8df$Stra8_Count_DEGs/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Unbound")$Stra8_Count_DEGs) * 100, 0))



#cgi degs
Stra8df$CGI_Count_DEGs <- c(nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Bound" & germ_KDM5C_Stra8_DEGs$Promoter_CpG_island == "noCGI")) , nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Bound" & germ_KDM5C_Stra8_DEGs$Promoter_CpG_island == "CGI")), nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8_DEGs$Promoter_CpG_island == "noCGI")), nrow(subset(germ_KDM5C_Stra8_DEGs, germ_KDM5C_Stra8_DEGs$KDM5C_binding == "Unbound" & germ_KDM5C_Stra8_DEGs$Promoter_CpG_island == "CGI")))

Stra8df$CGI_Percent_DEGs_plot <- ifelse(Stra8df$Kdm5c_binding == "Bound", Stra8df$CGI_Count_DEGs/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Bound")$CGI_Count_DEGs) * 100, ifelse(Stra8df$Kdm5c_binding == "Unbound", Stra8df$CGI_Count_DEGs/sum(subset(Stra8df, Stra8df$Kdm5c_binding == "Unbound")$CGI_Count_DEGs) * 100, 0))

Stra8df$Stra8_DEGs_Percent <- as.integer(round(Stra8df$Stra8_Percent_DEGs_plot))
Stra8df$CGI_DEGs_Percent <- as.integer(round(Stra8df$CGI_Percent_DEGs_plot))



Stra8df


#set the plotting order
motifbar <- function(colum, TITLE, FILL, COLORS){
	#df$Motifs <- factor(df$Motifs)

	
	library("ggpubr")
	q <- ggbarplot(Stra8df, "Kdm5c_binding", colum,
	fill = FILL, color = FILL, palette = COLORS,
	title = TITLE, label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "KDM5C Binding at Promoter", ylab = "% of genes", orientation = "vert") 

	return(q)
}

stra8colors <- c("Bound" = "firebrick1", "Unbound" = "brown4")

Stra8plots <- list()
Stra8plots[[1]] <- motifbar("Stra8_Percent", "All germline genes", "Stra8_binding", stra8colors)
Stra8plots[[2]] <- motifbar("Stra8_DEGs_Percent", "Germline DEGs", "Stra8_binding", stra8colors)

library("gridExtra")
ggsave(snakemake@output[[1]], grid.arrange(grobs = Stra8plots, ncol = 2), width = 8, height = 4)

cgicolors <- c("noCGI" = "aquamarine4", "CGI" = "aquamarine3")

CGI_plots <- list()
CGI_plots[[1]] <- motifbar("CGI_Percent", "All germline genes", "CpG_island", cgicolors)
CGI_plots[[2]] <- motifbar("CGI_DEGs_Percent", "Germline DEGs", "CpG_island", cgicolors)

library("gridExtra")
ggsave(snakemake@output[[2]], grid.arrange(grobs = CGI_plots, ncol = 2), width = 8, height = 4)
