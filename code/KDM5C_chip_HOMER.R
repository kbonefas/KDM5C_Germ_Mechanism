# 24.07.20 - plot the e2f and ebox motifs for kdm5c-bound and unbound germline genes

#read in everything
#make a dataframe that is bound and unbound together
#plot the data


#read in the HOMER output files for the histogram

motifs <- c("E2F", "E-box")
genes <- c("all germline genes", "germline DEGs")


#function to format the HOMER output
formatdf <- function(df){
	colnames(df) <- c("Distance_Center", "Total_Sites", "Pos_Sites", "Neg_Sites", "A_Freq", "C_Freq","G_Freq", "T_Freq")
	phist <- data.frame("TSS_Distance" = rep(df$Distance_Center, 2), "Sites_per_bp" = c(df$Pos_Sites, df$Neg_Sites), "Type" = c(rep("+ Sites", length(df$Pos_Sites)), rep("- Sites", length(df$Neg_Sites))))
	return(phist)
}


library('ggpubr')
histplot <- function(df, motif, n){
	p <- ggline(phist, "TSS_Distance", "Sites_per_bp", color = "Type", shape = 46, title = paste(motif, "binding sites")) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
	ggsave(snakemake@output[[n]], ggpar(p, legend = "bottom"), width = 6, height = 3.5)
}

#graphing bound and unbound on same graph
histplot_tog <- function(df, title, n){
	p <- ggline(df, "TSS_Distance", "Sites_per_bp", color = "Bound_Type", shape = 46, title =  title) + geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
	ggsave(snakemake@output[[n]], ggpar(p, legend = "right"), width = 9, height = 3.5)
}


makemotifs <- function(motif){
	offset <- ifelse(motif == "E2F", 0, ifelse(motif == "E-box", 4, "whoop"))

	for(i in 1:length(genes)){
		bound <- formatdf(read.csv(snakemake@input[[i+offset]], sep = "\t"))
		bound$binding <- "KDM5C bound"
		unbound <- formatdf(read.csv(snakemake@input[[i+length(genes)+offset]], sep = "\t"))
		unbound$binding <- "KDM5C unbound"

		tog <- rbind(bound,unbound)
		tog$Bound_Type <- paste(tog$binding, motif, tog$Type)
		print(head(tog))

		TITLE <- paste(motif, "binding sites at", genes[i])
		histplot_tog(tog, TITLE, i+(offset/length(genes)))
	}

}

for(m in motifs){
	makemotifs(m)
}