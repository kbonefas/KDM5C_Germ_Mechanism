#26.04.01 - Do kdm5c bound vs unbound correlate with stra8 targest

#Stra8 bound genes
Stra8 <- read.csv(snakemake@input[[2]], sep = ",")
Stra8 <- subset(Stra8, Stra8$STRA8_bound == "STRA8-bound")
#set the plotting order
motifbar <- function(colum, TITLE, FILL, COLORS){
	#df$Motifs <- factor(df$Motifs)

	
	library("ggpubr")
	q <- ggbarplot(Stra8df, "Kdm5c_binding", colum,
	fill = FILL, color = FILL, palette = COLORS,
	title = TITLE, label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "KDM5C Binding at Promoter", ylab = "% of genes", orientation = "vert") 

	return(q)
}

stra8colors <- c("yes" = "firebrick1", "no" = "brown4")
cgicolors <- c("no" = "aquamarine4", "CGI" = "aquamarine3")



############### stra8 and dazl, which are CGIs
#a and b are lists of genes
	#ex: a = stra8 targets
	#b = genes with CGIs
#labels = either "ENSEMBL" or "SYMBOL"

germ_percent_bar <- function(a, b, labels, NAME){
	#total number of genes
	together <- unique(c(a, b))
	print("together")
	print(head(together))

	a_in_b <- a[a %in% b]
	print("a_in_b")
	print(head(a_in_b))

	print("length a")
	print(length(a))

	print("length a_in_b")
	print(length(a_in_b))

	#instances of a that are in b
	#of all the germline genes that are stra8 targets (a), how many are in that have CGIs (b)
	a_in_b <- as.integer(round(length(a[a %in% b])/length(a)*100))
	print("")
	print(a_in_b)

	
	#how many that aren't stra8 targets have CGIs
	if(labels == "ENSEMBL"){
		IDs <- germ$ENSEMBL
		nota <- IDs[!(IDs %in% a)]
	} else if(labels == "SYMBOL"){
		IDs <- germ$SYMBOL
		nota <- IDs[!(IDs %in% a)]
	} else {
		print("whoops I broke")
	}
	print("nota")
	print(head(nota))

	nota_in_b <- as.integer(round(length(nota[nota %in% b])/length(nota)*100))

	### make the categories, name the categories, subtract percents
	df <- data.frame(CGI_status = c(rep("CGI", 2), rep("no", 2)),  Overlap = rep(c("yes", "no"), 2))
	df$percent <- c(a_in_b, 100-a_in_b, nota_in_b, 100-nota_in_b)

	print(df)


	
	q <- ggbarplot(df, "CGI_status", "percent",
	fill = "Overlap", color = "Overlap", palette = stra8colors,
	title = paste0(NAME, " overlap"), label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "CGI statusr", ylab = "% of genes", orientation = "vert") 

	return(q)
}

germ <- read.csv(snakemake@input[[1]], sep = ",")
CGI_genes <- subset(germ, germ$Promo_CGI == "CGI")


# (a, b, labels, NAME)


library("ggpubr")
ggsave(snakemake@output[[1]], germ_percent_bar(Stra8$ENSEMBL, CGI_genes$ENSEMBL, "ENSEMBL", "Stra8"), width = 4.5, height = 4)


