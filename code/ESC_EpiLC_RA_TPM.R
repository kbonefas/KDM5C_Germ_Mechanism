#24.03.15 Visualizing the TPM of genes as WT and Kdm5c cells differentiate from ESCs to EpiLCs with and without RA

#read in the tpm file
tpm_all <- read.csv(snakemake@input[[1]], sep ="\t", row.names = 1)
print(head(tpm_all[1:3]))
tpm <- tpm_all[4:ncol(tpm_all)] 
print("TPMs")
print(head(tpm))

#reorder the TPM  based on the samples
SampleInfo <- read.csv(snakemake@input[[2]], sep =",") 
rownames(SampleInfo) <- SampleInfo$ID

#set up the plotting variables
SampleInfo$RA[SampleInfo$RA == 'RA'] <- 'VA +'
SampleInfo$RA[SampleInfo$RA == 'NO'] <- 'VA -'

SampleInfo$genoRA <- paste(SampleInfo$Genotype, SampleInfo$RA)
SampleInfo$Timepoint <- as.factor(SampleInfo$Timepoint)

#make a line plot of the TPM of germline genes of interest
#geneID is the name of the gene symbol you want to plot

tpm_plot <- function(geneID){
		#get the ensembl ID for the gene you're plotting
	tpm_all$ENSEMBL <- rownames(tpm_all)
	ensembl <- tpm_all$ENSEMBL[tpm_all$external_gene_name == geneID]

	#make a dataframe with the sample timepoint, and genotype+RA
	germoi <- subset(tpm, rownames(tpm) == ensembl)
	TPMs <- t(germoi)
	# print(geneID)
	# print(head(TPMs))
	geneofinterest <- data.frame(TPM =TPMs[,1], ID = rownames(TPMs))
	geneofinterest <- merge(geneofinterest, SampleInfo, by = "ID")
	
	#EpiLCs with and without RA begin at the same starting point because they were differentiated from ESCs without RA 
		#make duplicate values for 0 hrs
	no0 <- subset(geneofinterest, RA != "")
	fix0 <- subset(geneofinterest, RA == "")
	fix0_RA <- fix0
	fix0_RA$genoRA <- paste(fix0_RA$Genotype, "VA +")
	fix0_NO <- fix0
	fix0_NO$genoRA <- paste(fix0_NO$Genotype, "VA -")
	fixed0 <- rbind(no0, fix0_RA, fix0_NO)

	#print(head(fixed0))
	source("code/utilities/colorpalettes.R")
	library("ggpubr")

	zerolab <- max(subset(fixed0, Timepoint == "0")$TPM)+0.5
	fortyeightlab <- max(subset(fixed0, Timepoint == "48")$TPM)+0.5
	ninetysixlab <-  max(subset(fixed0, Timepoint == "96")$TPM)+0.5

	#order the samples for plotting
	fixed0$genoRA <- factor(fixed0$genoRA , levels = c("WT VA +",  "5CKO VA +", "WT VA -",  "5CKO VA -"))

	my_comparisons <- c(list("5CKO VA +", "5CKO VA -"))
	#y coordinate for the significance star
	hjust <- 1.2 #adjustment factor
	y0 <- max(fixed0$TPM[fixed0$Timepoint == 0])*hjust
	y48 <- max(fixed0$TPM[fixed0$Timepoint == 48])*hjust
	y96 <- max(fixed0$TPM[fixed0$Timepoint == 96])*hjust
	#label.y = c(y0, y48, y96) 


	pgerm <- ggline(fixed0, x = "Timepoint", y = "TPM",  palette = genoRAcolors, shape = "genoRA", color = "genoRA", size = 1, point.size = 1.5, add = "mean_se", title = paste0(geneID)) 
	#stat_compare_means(comparisons = my_comparisons, label = "p.signif")

	pgerm <- ggpar(pgerm, ylim = c(0, ceiling(max(fixed0$TPM)*1.1)), legend = "right", font.main = "italic", xlab = "hrs differentiation", legend.title = "Genotype/RA")
	return(pgerm)

}



#just plotting 5CKO conditions
tpm_plot_5C <- function(geneID){
	#get the ensembl ID for the gene you're plotting
	tpm_all$ENSEMBL <- rownames(tpm_all)
	ensembl <- tpm_all$ENSEMBL[tpm_all$external_gene_name == geneID]

	#make a dataframe with the sample timepoint, and genotype+RA
	germoi <- subset(tpm, rownames(tpm) == ensembl)
	TPMs <- t(germoi)
	# print(geneID)
	# print(head(TPMs))
	geneofinterest <- data.frame(TPM =TPMs[,1], ID = rownames(TPMs))
	geneofinterest <- merge(geneofinterest, SampleInfo, by = "ID")
	
	#because plotting is dumb, add duplicate values for 0 hrs (one with RA and one without)
	no0 <- subset(geneofinterest, RA != "")
	fix0 <- subset(geneofinterest, RA == "")
	fix0_RA <- fix0
	fix0_RA$genoRA <- paste(fix0_RA$Genotype, "VA +")
	fix0_NO <- fix0
	fix0_NO$genoRA <- paste(fix0_NO$Genotype, "VA -")
	fixed0 <- rbind(no0, fix0_RA, fix0_NO)

	#print(head(fixed0))
	source("code/utilities/colorpalettes.R")
	library("ggpubr")

	zerolab <- max(subset(fixed0, Timepoint == "0")$TPM)+0.5
	fortyeightlab <- max(subset(fixed0, Timepoint == "48")$TPM)+0.5
	ninetysixlab <-  max(subset(fixed0, Timepoint == "96")$TPM)+0.5

	#order the samples for plotting
	fixed0$genoRA <- factor(fixed0$genoRA , levels = c("WT VA +",  "5CKO VA +", "WT VA -",  "5CKO VA -"))

	fixed0_5C <- subset(fixed0, fixed0$genoRA =="5CKO VA +" | fixed0$genoRA =="5CKO VA -")
	#y coordinate for the significance star
	hjust <- 1.05 #adjustment factor
	y0 <- max(fixed0_5C$TPM[fixed0_5C$Timepoint == 0])*hjust
	y48 <- max(fixed0_5C$TPM[fixed0_5C$Timepoint == 48])*hjust
	y96 <- max(fixed0_5C$TPM[fixed0_5C$Timepoint == 96])*hjust
	#label.y = c(y0, y48, y96)


	pgerm <- ggline(fixed0_5C, x = "Timepoint", y = "TPM",  palette = genoRAcolors, shape = "genoRA", color = "genoRA", size = 1, point.size = 1.5, add = "mean_se", title = paste0(geneID)) +
	stat_compare_means(aes(group = genoRA), method = "t.test", label = "p.signif", label.y = c(y0, y48, y96))

	pgerm <- ggpar(pgerm, ylim = c(0, ceiling(max(fixed0_5C$TPM)*1.1)), legend = "right", font.main = "italic", xlab = "hrs differentiation", legend.title = "Genotype/VA") 
	return(pgerm)

}







# NotRAgenies <- c("Dazl", "D1Pas1", "Sycp2", "Tex15", "Boll")
# RAgenies <- c("Stra8", "Ddx4", "Asz1", "Mael", "Sycp3", "Naa11", "Zcwpw1", "Mei1")

#geneielistie <- c("Ddx4", "Sycp3", "Tcfl5", "Stra8", "Mael", "Dazl", "Best1", "Ttc16", "Spag16", "D1Pas1")

#genes from each cluster
#2 - Pld6, Catsperg2
#3 - Dazl, D1pas1
#5 - Ddx4, Stra8
#7 - Spag16, Zpbp2


#1 - Tekt5, Dkkl1
#2 - Dazl, D1pas1
#3 - Boll, Tex16
#4 - Mael, Sycp3
#5 - Aste1, Tekt1
#6 - Spag16, Best1


geneielistie <- c("Pld6", "Catsperg2", "Dazl", "D1Pas1", "Sycp3", "Mael", "Spag16", "Zpbp2")

#empty lists to store the plots
exprTPM <- list()
exprTPM_5C <- list()

count <- 1
for (g in geneielistie){
	exprTPM[[count]] <- tpm_plot(g)
	exprTPM_5C[[count]] <- tpm_plot_5C(g)

	count <- count + 1
}

library("gridExtra")
print("making the tpm plots:")
ggsave(snakemake@output[[1]], plot = grid.arrange(grobs = exprTPM, nrow = 4), width = 8, height = 11)

ggsave(snakemake@output[[2]], plot = grid.arrange(grobs = exprTPM_5C, nrow = 4), width = 8, height = 11)
#warnings()



#genes of interest testing
geneielistie <- c("Mael", "Ddx4", "Piwil1", "Piwil2", "Asz1", "Tex15", "Tdrd1", "Tdrd9")

#empty lists to store the plots
exprTPM <- list()
count <- 1
for (g in geneielistie){
	exprTPM[[count]] <- tpm_plot(g)
	count <- count + 1
}

library("gridExtra")
print("making the tpm plots:")
ggsave(snakemake@output[[3]], plot = grid.arrange(grobs = exprTPM, nrow = 4), width = 8, height = 12)
#warnings()
