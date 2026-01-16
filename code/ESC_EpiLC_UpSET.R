#2023.01.13 - Upset plot germline genes 
germ <- read.csv(snakemake@input[[1]])
print("germline genes")
print(head(germ))

DEGs <- list()
#read in the DEGs
samples <- c("ESC", "48NO", "48VA", "96NO", "96VA")
for(i in 1:length(samples)){
    df <- read.csv(snakemake@input[[1+i]], sep = ",", row.names = 1)
    print(head(df))

    DEGs[[i]] <- row.names(df)
}

names(DEGs) <- samples


#### plot simple upset
library("UpSetR")

# modifiedupset <- function(samplelist){
# 	upset(fromList(samplelist), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "# of Overlapping Germline DEGs", empty.intersections = "on")
# }

pdf(file = snakemake@output[[1]], width = 8, height = 6)

upset(fromList(DEGs), order.by = "freq",  sets.x.label = "# germline DEGs", mainbar.y.label = "# in group", text.scale = 2, sets = samples, mb.ratio = c(0.55, 0.45), keep.order = TRUE)

dev.off()


#vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
#need a cagetory row in my upset list??? that has the category KDM5C bound
    #are you even able to do it in list format
# apparently you can do scatter plots with UpsetR
# Comparing two sets in a whole: https://github.com/hms-dbmi/UpSetR?tab=readme-ov-file



#plotting euler
library("ggplot2")



library("eulerr")
source("code/utilities/colorpalettes.R")


#subset the DEGs for just the RA treatment (aka "normal" differentiation)
    #just the unique ones?
DEGs_normal <- DEGs[c("ESC", "48VA", "96VA")] 
    #instead of doing Euler could do bar plot of the unique ones with darker shading for KDM5C binding - doesn't show overlap, would need to do another upset bar plot or something
    #take the unique DEGs and plot KDM5C binding vs shared 


KDM5C_bound <- subset(germ, germ$KDM5C_binding == "Bound")
DEGs_normal_KDM5C <- DEGs_normal
DEGs_normal_KDM5C[["KDM5C_bound"]] <- KDM5C_bound$ENSEMBL

together <- euler(DEGs_normal_KDM5C)

q <- plot(together, quantities = TRUE, labels = list(font = 4))
#fills = c(EpiLC_XY_KO, EpiLC_XX_HET, EpiLC_XX_KO)

library("ggplot2")
ggsave(snakemake@output[[2]], plot = q, width = 4, height = 4)

## only unique DEGs
#need unique DEGs
# DEGs_unique <- Reduce(setdiff, DEGs)
# print(DEGs_unique)

DEGs_unique <- list()


for(i in 1:length(samples)){
    #the samples that aren't the one you're testing
    others <- samples[samples != samples[i]]
    print(others)

    #get all the genes that are in the group not in your list
    others_colapse <- unlist(DEGs[others])
    
    DEGs_unique[[i]] <- setdiff(DEGs[[i]], others_colapse)

}

names(DEGs_unique) <- samples

print("DEGs_unique")

print(DEGs_unique)

#shared 
DEGs_unique[["shared"]] <- Reduce(intersect, DEGs)

# unique_96VA <- 
# 96only <- unique_within_list <- lapply(my_list_with_dups, unique)

DEGs_unique_kdm5c <- DEGs_unique

DEGs_unique_kdm5c[["KDM5C_bound"]] <- KDM5C_bound$ENSEMBL

together <- euler(DEGs_unique_kdm5c)

u <- plot(together, quantities = TRUE, labels = list(font = 4))
#fills = c(EpiLC_XY_KO, EpiLC_XX_HET, EpiLC_XX_KO)

ggsave(snakemake@output[[3]], plot = u, width = 4, height = 4)


###upset plot of normal genes
pdf(file = snakemake@output[[4]], width = 8, height = 6)

upset(fromList(DEGs_unique), order.by = "freq",  sets.x.label = "# germline DEGs", mainbar.y.label = "# in group", text.scale = 2, sets = c("ESC", "48VA", "96VA"), mb.ratio = c(0.55, 0.45), keep.order = TRUE)

dev.off()




################### germ CGI stage 
#expression of genes in spermatogenesis stages
greengerm <- read.csv(snakemake@input[[7]], sep = ",")

#plotting df columns:
	#germ cell stage
	#transient or silent
	#average expression (for all transient or silent)
green_plot <- data.frame()
expression <- c(samples, "shared")
print(expression)

for(i in expression){
	germ_gene <- subset(germ, germ$ENSEMBL %in% DEGs_unique[[i]])

    symbols <- germ_gene$SYMBOL
	#for every gene symbol, get the expression value and add it to the dataframe
	for(k in symbols){
		expr <- subset(greengerm, greengerm$SYMBOL == k)
	
		germexpr <- c(t(expr[1,2:ncol(expr)])) #values of the expression
		
		#get the name of the germ cell stage and the means of expression, skipping the first column (gene symbols)
		df <- data.frame(Stage = colnames(greengerm)[2:ncol(greengerm)], GermExpr = germexpr, Group = rep(i, ncol(greengerm)-1))

		green_plot <- rbind(green_plot, df)

	}
		
}

print(head(green_plot))

library('ggpubr')

my_comparisons <- expression
p <- ggboxplot(green_plot, "Stage", "GermExpr", fill = "Group",  title = "Average expression in germ cell stages", ylab = "log(Avg of Normalized Expression + 1)", xlab = "Stage of Spermatogenesis", outlier.shape = NA, ylim = c(0,5)) +
	stat_compare_means(aes(label=..p.signif.., group=Group), method="wilcox.test", label.y = 4, size = 8) +
	font("xy.text", size = 24) + font("title", size = 35, face = "bold") + font("xlab", size = 24) + font("ylab", size = 24)

#palette = c("#f93a0b", "#ff8a7a"),
p <- ggpar(p, legend = "top", legend.title = "DEG group")

ggsave(snakemake@output[[5]], p, width = 20, height = 8)


## heatmap of expression across stages
	#multiple heatmaps of each developmental stage stacked on top 
#make a hclust matrix for each stage, have the genes expression in CGI vs no. 

#for every germ cell stage, get the germ cell expression and for CGI and CGI no genes

#make 2 heatmaps (CGI vs no), split by germ cell stage. Use the same scale for both so they can be compared to each other

#each column is a germ cell stage, each row is a gene, values are the expression
library(ComplexHeatmap)
library(circlize)

expression <- c("CGI", "no")
#setting the range of colors
heat_color = colorRamp2(range(subset(greengerm, select = c(-SYMBOL))), hcl_palette = "Reds", reverse = TRUE)


for(i in 1:length(expression)){
	symbols <- subset(germ, Promo_CGI == expression[i])[,"SYMBOL"]

	df <- subset(greengerm, greengerm$SYMBOL %in% symbols)
	row.names(df) <- df$SYMBOL
	plotmatrix <- subset(df, select = c(-SYMBOL))
	p <- Heatmap(plotmatrix, show_row_names = FALSE, column_title = paste("Germline genes with", expression[i]), cluster_columns = FALSE, heatmap_legend_param = list(title = "log(Avg of Normalized Expression + 1)"), col = heat_color)
	

	pdf(file = snakemake@output[[i+1]],   # The directory you want to save the file in
    	width = 10, # The width of the plot in inches
    	height = 10) # The height of the plot in inches
		draw(p)
	dev.off()
		
}

















############### stacked barplot of kdm5c binding unique DEGs
# allgerm <- subset(plotdf, plotdf$Kdm5c_binding == "All germ")
# allgerm <- subset(allgerm, select = c("Kdm5c_binding", "CpG_island", "CGIPercent_plot", "Percent"))
# colnames(allgerm) <- c("Gene_type", "CpG_island", "raw", "Percent")
# plotdf2 <- rbind(plotdf2, allgerm)

#columns - time point, KDM5C binding status (KDM5C_binding), #, percent
# one column - 96hrs, Kdm5c_bound, # bound, % of total
# 2nd column - 96hrs, Kdm5c_unbound, # unbound, %  of total




#plot the results in a bar graph
#set plotting order
# plotdf2$Gene_type <- factor(plotdf2$Gene_type, levels = c("All genes", "All germ"))

# library("ggpubr")
# all_p <- ggbarplot(plotdf2, "Gene_type", "Percent", fill = "CpG_island", color = "CpG_island", palette = c("no" = "#ff8a7a", "CGI" = "#f93a0b"),
# 		title = "CpG islands at mm10 promoters", label = TRUE, lab.col = "white", lab.vjust = 1, xlab = " ", ylab = "% of genes", orientation = "vert") 

# ggsave(snakemake@output[[8]], all_p, width = 3, height = 4)















# library(ComplexUpset)
# movies = as.data.frame(ggplot2movies::movies)
# genres = colnames(movies)[18:24]

# # for simplicity of examples, only use the complete data points
# movies[movies$mpaa == '', 'mpaa'] = NA
# movies = na.omit(movies)


# upset(
#     movies,
#     genres,
#     base_annotations=list(
#         'Intersection size'=intersection_size(
#             counts=FALSE,
#             mapping=aes(fill=mpaa)
#         )
#     ),
#     width_ratio=0.1
# )

# library(UpSetR)
# test <- as.data.frame(matrix(rnorm(1:12), nrow = 4))
# colnames(test) <- c("A", "B", "C")
# test_up <- test
# test_up[test >= 0] <- 1
# test_up[test_up != 1] <- 0
# test_up$Gene <- paste0("Gene", 1:nrow(test), "_UP")
# test_down <- test
# test_down[test < 0] <- 1
# test_down[test_down != 1] <- 0
# test_down$Gene <- paste0("Gene", 1:nrow(test), "_DOWN")
# test <- rbind(test_up, test_down)


# upordown <- function(row, direction) {
#   gene <- row["Gene"]
#   if (grepl(x = gene, pattern = direction)) {
#     newData <- T
#   }
#   else {
#     newData <- F
#   }
# }


# metadata <- data.frame(
#   c("A", "B", "C"),
#   as.numeric(apply(test_up[, 1:3], 2, sum))
# )
# colnames(metadata) <- c(
#   "sets",
#   "NumberUP"
# )


# upset(test,
#       sets = c("A", "B", "C"), set.metadata = list(
#         data = metadata,
#         plots = list(
#           list(type = "hist", 
#                column = "NumberUP", 
#                assign = 20, # defines width of the meta-data histogram
#                colors = "red")
#         )
#       ),
#       queries = list(list(query = upordown, params = list("_UP"), color = "red", active = TRUE))
# )