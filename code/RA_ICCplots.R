#2025.03.18 - Plots for ICC counts in WT and 5CKO cells
library("readxl")
DAZL_100nMRA <- read_excel(snakemake@input[[1]])

#make a boxplot
library("ggplot2")
library("ggpubr")

#make a new column with genotype and treatment to compare
DAZL_100nMRA$genotreat <- paste(DAZL_100nMRA$Genotype, DAZL_100nMRA$Treatment)
#make a new column with the percent of DAPI cells that are DAZL+
DAZL_100nMRA$Percent <- (DAZL_100nMRA$DAZL/DAZL_100nMRA$DAPI*100)

print(head(DAZL_100nMRA))




#Plot data based on average of each sample
#need to include the collection date and split by RA and DMSO treatment

#list of sample names
DAZL_100nMRA$ID <- paste(DAZL_100nMRA$Sample, DAZL_100nMRA$genotreat)
allsamps <- unique(DAZL_100nMRA$ID)
dflist <- list()

for (k in 1:length(allsamps)){
    #get the dataframe for one sample
    df <- subset(DAZL_100nMRA, DAZL_100nMRA$ID == allsamps[k])
    #print(df)
    newdf <- df[1,3:ncol(df)]
    print(newdf)

    #newdf$ID <- allsamps[k]
    newdf$Avg_perc <- mean(df$Percent)
	dflist[[k]] <- newdf

}

library(dplyr)
DAZL_avg <- bind_rows(dflist)

print(head(DAZL_avg))

#round to the nearest 2 decimal places
#round(DATA, digits = 2)


##### Plotting
#order the samples for plotting
DAZL_100nMRA$genotreat <- factor(DAZL_100nMRA$genotreat, levels = c("WT DMSO", "5CKO DMSO", "WT RA", "5CKO RA"))


my_comparisons <- list(c("WT RA", "5CKO RA"), c("WT DMSO", "5CKO DMSO"), c("5CKO DMSO", "5CKO RA"), c("WT DMSO", "WT RA"))

source("code/utilities/colorpalettes.R")
q <- ggbarplot(DAZL_avg, x = 'genotreat', y = 'Avg_perc', fill="genotreat", add = c("mean_se", "jitter"),
    xlab = " ", ylab = "% DAZL+/DAPI+", palette = genoRAcolors) + 
    #rremove("legend") +
    stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif") 
q <- ggpar(q, x.text.angle = 25, font.main = "bold", legend = "right", legend.title = " ")

# q <- facet(q, facet.by = "Symbol", nrow = 2)

ggsave(snakemake@output[[1]], plot = q, width = 5, height = 4)

