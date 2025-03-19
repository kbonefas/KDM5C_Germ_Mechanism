#plotting number of tissue-enriched genes are in each Kdm5c-KO dataset 


#1) Get the number and p-value of tissue-enriched DEGs for each Kdm5c-KO sample 
#samples <- c("nESC", "EpiLC", "exEpiLC", "Amygdala", "Hippocampus")
samples <- c("Amygdala", "Hippocampus")

source('code/utilities/tissueSPECIFICgenes.R')
print(head(Tissues))


tissuecount <- data.frame()

dflist <- list()
for(i in 1:length(samples)){
	dflist[[i]] <- tissue_dot(i, samples[i])
}

library(dplyr)
tissuecount <- bind_rows(dflist)

print(head(tissuecount))

#2) set samples with zero tissue DEGs to NA so they plot as blank
tissuecount$DEGnumb[tissuecount$DEGnumb == 0] <- NA

#get lowest p-value for the minimum limit
mingrad <- min(tissuecount$pvalue)
print(mingrad)

#3) Plot
library("ggplot2")
p <- ggplot(tissuecount, aes(x=tissue, y = sample, color = pvalue, size = DEGnumb)) + labs(x = NULL, y = NULL) +
  geom_point() + scale_color_gradient(low="red", high="blue", limits = c(mingrad, 0.05)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(snakemake@output[[1]], plot = p, width = 7, height = 4)
