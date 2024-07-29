# 24.07.20 - plot the e2f and ebox motifs for kdm5c-bound and unbound germline genes

#read in everything
#make a dataframe that is bound and unbound together
#plot the data


# #read in the HOMER output files for the histogram

# motifs <- c("E2F", "E-box")
# genes <- c("all germline genes", "germline DEGs")


# #function to format the HOMER output
# formatdf <- function(df){
# 	colnames(df) <- c("Distance_Center", "Total_Sites", "Pos_Sites", "Neg_Sites", "A_Freq", "C_Freq","G_Freq", "T_Freq")
# 	phist <- data.frame("TSS_Distance" = rep(df$Distance_Center, 2), "Sites_per_bp" = c(df$Pos_Sites, df$Neg_Sites), "Type" = c(rep("+ Sites", length(df$Pos_Sites)), rep("- Sites", length(df$Neg_Sites))))
# 	return(phist)
# }


# library('ggpubr')
# histplot <- function(df, motif, n){
# 	p <- ggline(phist, "TSS_Distance", "Sites_per_bp", color = "Type", shape = 46, title = paste(motif, "binding sites")) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
# 	ggsave(snakemake@output[[n]], ggpar(p, legend = "bottom"), width = 6, height = 3.5)
# }

# #graphing bound and unbound on same graph
# histplot_tog <- function(df, title, n){
# 	p <- ggline(df, "TSS_Distance", "Sites_per_bp", color = "Bound_Type", shape = 46, title =  title) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
# 	ggsave(snakemake@output[[n]], ggpar(p, legend = "right"), width = 9, height = 3.5)
# }


# makemotifs <- function(motif){
# 	offset <- ifelse(motif == "E2F", 0, ifelse(motif == "E-box", 4, "whoop"))

# 	for(i in 1:length(genes)){
# 		bound <- formatdf(read.csv(snakemake@input[[i+offset]], sep = "\t"))
# 		bound$binding <- "KDM5C bound"
# 		unbound <- formatdf(read.csv(snakemake@input[[i+length(genes)+offset]], sep = "\t"))
# 		unbound$binding <- "KDM5C unbound"

# 		tog <- rbind(bound,unbound)
# 		tog$Bound_Type <- paste(tog$binding, motif, tog$Type)
# 		print(head(tog))

# 		TITLE <- paste(motif, "binding sites at", genes[i])
# 		histplot_tog(tog, TITLE, i+(offset/length(genes)))
# 	}

# }

# for(m in motifs){
# 	makemotifs(m)
# }

# 24.07.28 - percentage of promoters with e2f and ebox  motifs

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

motifdf <- data.frame(Kdm5c_binding = c(rep("Bound", 4), rep("Unbound", 4)), Motifs = c("E2F", "Ebox", "Both", "Neither"))
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


motifdf$Percent <- ifelse(motifdf$Kdm5c_binding == "Bound", (motifdf$Count/(length(E2F_only_bound) + length(Ebox_only_bound) + length(Both_bound) + length(Neither_bound))) * 100, ifelse(motifdf$Kdm5c_binding == "Unbound", (motifdf$Count/(length(E2F_only_unbound) + length(Ebox_only_unbound) + length(Both_unbound) + length(Neither_unbound))) * 100, "oops"))


motifdf
#save the plots in an empty list
plots <- list()

#set the plotting order
motifbar <- function(df, colum, TITLE){
	df$Motifs <- factor(df$Motifs, levels = c("E2F", "Ebox", "Both", "Neither"))
	library("ggpubr")
	q <- ggbarplot(df, "Kdm5c_binding", colum,
	fill = "Motifs", color = "Motifs", palette = c("E2F" = "forestgreen", "Ebox" = "blue3", "Both" = "goldenrod2", "Neither" = "slategray4"),
	title = TITLE, label = TRUE, lab.col = "black", lab.vjust = 1, xlab = "RNA-seq germline DEG", ylab = "# of germline DEGs", orientation = "vert") 

	return(q)
}

plots[[1]] <- motifbar(motifdf, "Count", "All germline genes")
plots[[2]] <- motifbar(motifdf, "Percent", "All germline genes")


library("gridExtra")
ggsave(snakemake@output[[1]], grid.arrange(grobs = plots, ncol = 2), width = 6.5, height = 3.5)



# e2finst <- read.csv(snakemake@input[[9]], sep = "\t")
# print(head(e2finst))

