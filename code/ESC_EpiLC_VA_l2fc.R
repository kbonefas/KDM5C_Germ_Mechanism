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

##Stra8 bound genes
Stra8 <- read.csv(snakemake@input[[6]], sep = ",")

#DAZL regulated genes
library("readxl")
Dazl <- data.frame(read_excel(snakemake@input[[7]], sheet = 5))
Dazl <- subset(Dazl, Dazl$DAZL.target == "DAZL.target")
print("Dazl")
print(head(Dazl))


l2fc_scatter <- function(df, TITLE){
    colnames(df)[colnames(df) == 'Row.names'] <- 'ENSEMBL'
    df2 <- merge(df, germ, by = "ENSEMBL") 
    df2$Stra8 <- ifelse(df2$ENSEMBL %in% Stra8$ENSEMBL, "yes", "no" )
    df2$Dazl <- ifelse(df2$SYMBOL %in% Dazl$Gene.id, "Dazl", "no" )
    print(head(df2))
    p <- ggscatter(df2, x = "NO_L2FC", y = "VA_L2FC",
        color = "Dazl", palette = c("blue","gray"), size = 2.5, alpha = 0.25, # Points color, shape and size
        add = "reg.line",  # Add regressin line
        add.params = list(color = "red", fill = "gray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
        title = TITLE
        ) + geom_abline(intercept = 0, slope = 1, linetype="solid", size=0.75)

    p + geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
}

ggsave(snakemake@output[[1]], ggpar(l2fc_scatter(hrs48, "48hrs"), xlim = c(-2,4), ylim = c(-2,4)), width = 5, height = 5)

ggsave(snakemake@output[[2]], ggpar(l2fc_scatter(hrs96, "96hrs"), xlim = c(-5,10), ylim = c(-5,10)), width = 5, height = 5)

