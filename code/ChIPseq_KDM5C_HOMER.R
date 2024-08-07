# 24.07.28 - percentage of KDM5C-bound and unbound promoters with e2f and ebox motifs
library("ggpubr")

#read in the KDM5C bound genes
germ_KDM5C <- read.csv(snakemake@input[[1]], sep = ",")
print(head(germ_KDM5C))
#column of kdm5c binding
#KDM5C_binding

#info we need
#total number bounds vs unbound kdm5c
#number with E2F
#number with ebox
#number with both
#number with neither

#make a column that is E2F and the answer is y/n. Then make a final colum that's motif status to compile the info

#plotting df - 
	#Kdm5c_binding
	#Motifs
	#Percentage

motifdf <- data.frame(Kdm5c_binding = c(rep("Bound", 4), rep("Unbound", 4)), Motifs = c("E2F", "E-box", "Both", "Neither"))
print(motifdf)


#read in the instances

# get unique ensembl ids

#e2f bound
E2F_bound <- unique(read.csv(snakemake@input[[2]], sep = "\t")$Ensembl)
#ebox bound
Ebox_bound <- unique(read.csv(snakemake@input[[3]], sep = "\t")$Ensembl)

Both_bound <- E2F_bound[E2F_bound %in% Ebox_bound]
E2F_only_bound <- E2F_bound[!E2F_bound %in% Ebox_bound]
Ebox_only_bound <- Ebox_bound[!Ebox_bound %in% E2F_bound]

germ_KDM5C_bound <- subset(germ_KDM5C, germ_KDM5C$KDM5C_binding == "Bound")$ENSEMBL
Neither_bound <- germ_KDM5C_bound[!(germ_KDM5C_bound %in% E2F_bound & germ_KDM5C_bound %in% Ebox_bound)]



#e2f unbound
E2F_unbound <- unique(read.csv(snakemake@input[[4]], sep = "\t")$Ensembl)
#ebox unbound
Ebox_unbound <- unique(read.csv(snakemake@input[[5]], sep = "\t")$Ensembl)

Both_unbound <- E2F_unbound[E2F_unbound %in% Ebox_unbound]
E2F_only_unbound <- E2F_unbound[!E2F_unbound %in% Ebox_unbound]
Ebox_only_unbound <- Ebox_unbound[!Ebox_unbound %in% E2F_unbound]

germ_KDM5C_unbound <- subset(germ_KDM5C, germ_KDM5C$KDM5C_binding == "Unbound")$ENSEMBL
Neither_unbound <- germ_KDM5C_unbound[!germ_KDM5C_unbound %in% E2F_unbound & !germ_KDM5C_unbound %in% Ebox_unbound]


motifdf$Count <- c(length(E2F_only_bound), length(Ebox_only_bound), length(Both_bound), length(Neither_bound), length(E2F_only_unbound), length(Ebox_only_unbound), length(Both_unbound), length(Neither_unbound))


motifdf$Percent_plot <- ifelse(motifdf$Kdm5c_binding == "Bound", (motifdf$Count/(length(E2F_only_bound) + length(Ebox_only_bound) + length(Both_bound) + length(Neither_bound))) * 100, ifelse(motifdf$Kdm5c_binding == "Unbound", (motifdf$Count/(length(E2F_only_unbound) + length(Ebox_only_unbound) + length(Both_unbound) + length(Neither_unbound))) * 100, 0))


motifdf$Percent <-as.integer(round(motifdf$Percent_plot))


#save the plots in an empty list
plots <- list()

library(wesanderson)

#set the plotting order

motifbar <- function(df, colum, TITLE){
	df$Motifs <- factor(df$Motifs, levels = c("E2F", "E-box", "Both", "Neither"))
	q <- ggbarplot(df, "Kdm5c_binding", colum,
	fill = "Motifs", color = "Motifs", palette = wes_palette("GrandBudapest2", 4),
	title = TITLE, label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "KDM5C Binding at Promoter", ylab = "% of genes with motif", orientation = "vert") 

	return(q)
}


ggsave(snakemake@output[[1]], motifbar(motifdf, "Percent", "All germline genes"), width = 4, height = 4)


#X-box motifs
xbox_df <- data.frame(Kdm5c_binding = rep(c("Bound", "Unbound"), 2), Xbox_status = c(rep("X-box",2), rep("No", 2)))

#####
xbox_KDM5C_bound <- unique(read.csv(snakemake@input[[6]], sep = "\t")$Ensembl)
xbox_KDM5C_unbound <- unique(read.csv(snakemake@input[[7]], sep = "\t")$Ensembl)

No_KDM5C_bound <- germ_KDM5C_bound[!germ_KDM5C_bound %in% xbox_KDM5C_bound]
No_KDM5C_unbound <- germ_KDM5C_unbound[!germ_KDM5C_unbound %in% xbox_KDM5C_unbound]

xbox_df$Count <- c(length(xbox_KDM5C_bound), length(xbox_KDM5C_unbound), length(No_KDM5C_bound), length(No_KDM5C_unbound))

xbox_df$Percent_plot <- ifelse(xbox_df$Kdm5c_binding == "Bound", xbox_df$Count/sum(subset(xbox_df, xbox_df$Kdm5c_binding == "Bound")$Count) * 100, ifelse(xbox_df$Kdm5c_binding == "Unbound", xbox_df$Count/sum(subset(xbox_df, xbox_df$Kdm5c_binding == "Unbound")$Count) * 100, 0))

xbox_df$Percent <- as.integer(round(xbox_df$Percent_plot))

xbox_df

xbox_df$Xbox_status <- factor(xbox_df$Xbox_status, level = c("X-box", "No"))

q <- ggbarplot(xbox_df, "Kdm5c_binding", "Percent", fill = "Xbox_status", color = "Xbox_status", palette = c("royalblue1", "royalblue4"),
title = "All germline genes", label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "KDM5C Binding at Promoter", ylab = "% of genes", orientation = "vert") 

ggsave(snakemake@output[[2]], q, width = 4, height = 4)



