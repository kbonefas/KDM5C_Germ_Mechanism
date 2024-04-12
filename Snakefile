
########################################################
#                                                      #
#     Snakefile for KDM5C Germline Mechanism Paper     #
#                                                      #
########################################################

####################### Figure 1 #######################

##padj cutoff for DESEq2
PADJ = 0.1
##number of genotypes you're comparing to WT

#generate DESeq2 results for 5CKO amygdala and hippocampus
rule amyhipDESeq2:
	input:
		#counts file
		"data/raw/20221005_adultHIPAMY_Kdm5cKO_cts.txt",
		"data/raw/SampleInfo_amyhip.csv" 
	params:
		alpha = PADJ, #padj cutoff
	output:
		"results/figure_pieces/PCA_HIPAMY.pdf", #PCA plot
		#the results tables will be generated in alphabetical order 
		"data/processed/restable_amy5cKO.csv",
		"data/processed/restable_hip5cKO.csv",
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv",
		"results/figure_pieces/GOmap_amy.pdf",
		"results/figure_pieces/GOmap_hip.pdf",
		"results/figure_pieces/GOdot_amy.pdf",
		"results/figure_pieces/GOdot_hip.pdf"
	script:
		"code/DESeq2_adultHIPAMY.R"


#look at tissue-enriched gene expression
rule tissue_genes:
	input:
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv",
		"data/processed/restable_amy5cKO.csv",
		"data/processed/restable_hip5cKO.csv"

	output:
		"results/figure_pieces/TissueSpecific_Volcano_amy5cKO.png",
		"results/figure_pieces/TissueSpecific_Volcano_hip5cKO.png",
		"results/figure_pieces/TissueSpecific_bar_amy5cKO.pdf",
		"results/figure_pieces/TissueSpecific_bar_hip5cKO.pdf",
		"results/figure_pieces/TissueSpecific_MAplot_amy5cKO.pdf",
		"results/figure_pieces/TissueSpecific_MAplot_hip5cKO.pdf",
		"data/processed/TissueSpecific_AMYHIP_numberofgenesDEGs.csv",
		"data/processed/TestisDEGs_amyhip.csv",
		"results/figure_pieces/TissueSpecific_bar_numbgenes.pdf"
	script:
		"code/TissueGenes_AmyHip.R"
