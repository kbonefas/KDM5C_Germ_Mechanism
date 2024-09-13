#24.09.13 Germline stages of development for CGI genes

#list of germline genes
germ <- read.csv(snakemake@input[[1]], sep = ",")


#expression of genes in spermatogenesis stages
greengerm <- read.csv(snakemake@input[[2]], sep = ",")



#plotting df columns:
	#germ cell stage
	#transient or silent
	#average expression (for all transient or silent)
green_plot <- data.frame()
expression <- c("CGI", "no")
for(i in expression){
	symbols <- subset(germ, Promo_CGI == i)[,"SYMBOL"]

	#for every symbol, get the expression value and add it to the dataframe
	for(k in symbols){
		expr <- subset(greengerm, greengerm$SYMBOL == k)
	
		germexpr <- c(t(expr[1,2:ncol(expr)])) #values of the expression
		
		#get the name of the germ cell stage and the means of expression, skipping the first column (gene symbols)
		df <- data.frame(Stage = colnames(greengerm)[2:ncol(greengerm)], GermExpr = germexpr, CGI_status = rep(i, ncol(greengerm)-1))

		green_plot <- rbind(green_plot, df)

	}
		
}

print(head(green_plot))

library('ggpubr')

my_comparisons <- expression
p <- ggboxplot(green_plot, "Stage", "GermExpr", fill = "CGI_status",  palette = c("#f93a0b", "#ff8a7a"), title = "Average expression in germ cell stages", ylab = "log(Avg of Normalized Expression + 1)", xlab = "Stage of Sperm Development") +
	stat_compare_means(aes(label=..p.signif.., group=CGI_status), method="wilcox.test", label.y = 6.5, size = 8) +
	font("xy.text", size = 24) + font("title", size = 35, face = "bold") + font("xlab", size = 24) + font("ylab", size = 24)

p <- ggpar(p, legend = "top", legend.title = "CGI at promoter")

ggsave(snakemake@output[[1]], p, width = 15, height = 10)
