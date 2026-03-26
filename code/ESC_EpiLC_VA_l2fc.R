#plot l2fc with RA vs no RA

#read in germ
germ <- read.csv(snakemake@input[[1]], sep = ",")

#make a dataframe 
times <- c("48 hrs", "96 hrs")

VA48 <- read.csv(snakemake@input[[2]], row.names = 1)
print(head(VA48))
VA48$VA_L2FC <- VA48$log2FoldChange
NO48 <- read.csv(snakemake@input[[3]], row.names = 1)
NO48$NO_L2FC <- NO48$log2FoldChange

VA96 <- read.csv(snakemake@input[[4]], row.names = 1)
VA96$VA_L2FC <- VA96$log2FoldChange
NO96 <- read.csv(snakemake@input[[5]], row.names = 1)
NO96$NO_L2FC <- NO96$log2FoldChange

hrs48 <- merge(VA48, NO48, by = 'row.names')
hrs96 <- merge(VA96, NO96, by = 'row.names')

print(head(hrs48))

#make a scatterplot
library(ggpubr)

l2fc_scatter <- function(df, TITLE){
    colnames(df)[colnames(df) == 'Row.names'] <- 'ENSEMBL'
    df2 <- merge(df, germ, by = "ENSEMBL") 
    print(head(df2))
    ggscatter(df2, x = "NO_L2FC", y = "VA_L2FC",
        color = "black", size = 2.5, # Points color, shape and size
        add = "reg.line",  # Add regressin line
        add.params = list(color = "red", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
        title = TITLE
        ) + geom_abline(intercept = 0, slope = 1, linetype="dashed", size=0.75)
}

ggsave(snakemake@output[[1]], ggpar(l2fc_scatter(hrs48, "48hrs"), xlim = c(-2,4), ylim = c(-2,4)), width = 5, height = 5)

ggsave(snakemake@output[[2]], ggpar(l2fc_scatter(hrs96, "96hrs"), xlim = c(-5,10), ylim = c(-5,10)), width = 5, height = 5)

