#2025.03.18 - Plots for ICC counts in WT and 5CKO cells
library("readxl")
DAZL_100nMRA <- read_excel(snakemake@input[[1]])

#make a boxplot
library("ggplot2")
library("ggpubr")

#make a new column with genotype and treatment to compare
DAZL_100nMRA$genotreat <- paste(DAZL_100nMRA$Genotype, DAZL_100nMRA$Treatment)
#make a new column with the percent of DAPI cells that are DAZL+
    #round to the nearest 2 decimal places
DAZL_100nMRA$Percent <- round(DAZL_100nMRA$DAZL/DAZL_100nMRA$DAPI*100, digits = 2)

print(head(DAZL_100nMRA))

my_comparisons <- list(c("WT RA", "5CKO RA"), c("WT DMSO", "5CKO DMSO"), c("5CKO DMSO", "5CKO RA"), c("WT DMSO", "WT RA"))

q <- ggbarplot(DAZL_100nMRA, x = 'genotreat', y = 'Percent', fill="genotreat", color = "genotreat", add = "mean_se",
    xlab = " ", ylab = "% DAZL+/DAPI+") + #palette = EpiLC_XY_palette
    #rremove("legend") +
    stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif") 
q <- ggpar(q, x.text.angle = 25, font.main = "bold")

# q <- facet(q, facet.by = "Symbol", nrow = 2)

ggsave(snakemake@output[[1]], plot = q, width = 6, height = 5)

