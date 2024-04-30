
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
rule brain_tissue_genes:
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
		"results/figure_pieces/TissueSpecific_bar_numbgenes.pdf",
		"data/processed/TissueDEGs_amy5CKO.csv",
		"data/processed/TissueDEGs_hip5CKO.csv"
	script:
		"code/TissueGenes_AmyHip.R"


###################### Figure 2: Germline-enriched Genes ############################
#gene ontology analysis of testis-enriched DEGs
rule testisGO:
	input:
		"data/processed/TestisDEGs_amyhip.csv"
	output:
		"results/figure_pieces/TestisDEGs_GO.pdf",
		"results/GO_Brain_TestisDEGs.csv"
	script:
		"code/GO_brain_testisDEGs.R"

#test if DEGs are specific to germ cells or somatic cells in the testis
rule testisEXPR:
	input:
		"data/raw/Mueller2013_adultWW_FPKM.txt", #FPKM adult WWv and WT testis
			#raw data from https://pubmed.ncbi.nlm.nih.gov/23872635/  
		"data/raw/Green2018_Testis_scRNAseq_Cell_Types.txt",
			#raw data from https://pubmed.ncbi.nlm.nih.gov/30146481/
		"data/processed/TestisDEGs_amyhip.csv"
	output:
		"results/figure_pieces/testisEXPR_mueller.pdf",
		"results/figure_pieces/testisEXPR_green.pdf"
	script:
		"code/testisEXPR_mueller_green.R"


#plotting ovary DEGs in germline-depleted ovary 
#rule ovaryEXPR:


#generate list of germline-enriched genes based on expression cutoff
rule germgenes:
	input:
		"data/raw/Mueller_embryonicmouseWWvsequencing.xlsx", #FPKM of embryonic WWv (germline-depeleted) and WT testis
		"data/raw/Mueller2013_adultWW_FPKM.txt", #FPKM adult WWv and WT testis
			#raw data from https://pubmed.ncbi.nlm.nih.gov/23872635/  
		"data/raw/Li_2017_Tissues_FPKMs.xlsx" #WT mouse tissue FPKM  
			#raw data from https://pubmed.ncbi.nlm.nih.gov/28646208/ 
	output:
		"data/processed/germGENES20.csv", #high confidence germline genes that lose 80% of their maximum expression in WT somatic tissues and germline-depleted testis
		#make graph of the total number of genes filtered to get the germline genes
		"results/figure_pieces/germGENESfilter.pdf"
	script:
		"code/utilities/mouseGermlineGenes.R"


###################### Figure 3: EpiLC germline genes ############################
#DESeq2 on Kdm5c-KO EpiLCs
rule EpiDESeq2:
	input:
		cts = "data/raw/EpiLC_gene_expected_count.txt"
	output:
		"data/raw/SampleInfo_EpiLC.csv",
		"data/processed/restable_EpiLC_XY5cKO.csv",
		"results/DESeq2/DEGs_EpiLC_XY5cKO.csv"
	params:
		alpha = PADJ
	script:
		"code/DESeq2_EpiLC.R"

#plot the expression of ESC and EpiLC pluripotency markers in EpiLCs
rule EpiLC_markers:
	input:
		"data/raw/EpiLC_gene_TPM.txt",
		"data/raw/SampleInfo_EpiLC.csv"
	output:
		"results/figure_pieces/EpiLC_markers_box_Supplement.pdf",
		"results/figure_pieces/EpiLC_markers_box.pdf"
	script:
		"code/EpiLC_markers.R"

#expression of tissue-specific genes in EpiLCs
rule EpiLC_tissue_genes:
	input:
		"results/DESeq2/DEGs_EpiLC_XY5cKO.csv",
		"data/processed/restable_EpiLC_XY5cKO.csv"

	output:
		"results/figure_pieces/TissueSpecific_Volcano_EpiLC_XY5cKO.png",
		"results/figure_pieces/TissueSpecific_bar_EpiLC_XY5cKO.pdf",
		"results/figure_pieces/TissueSpecific_MAplot_EpiLC_XY5cKO.pdf",
		"data/processed/TissueSpecific_EpiLC_XY_numberofgenesDEGs.csv",
		"data/processed/TissueDEGs_EpiLC5CKO_XY.csv"
	script:
		"code/TissueGenes_EpiLC.R"


#pull out germline DEGs in Kdm5c-KO EpiLCs and brain
rule BrainEpi_germDEGs:
	input:
		"data/processed/germGENES20.csv",
		"results/DESeq2/DEGs_EpiLC_XY5cKO.csv",
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv"
	output: #make sure match the order of input files
		"results/DESeq2/germDEGs/germDEGs_EpiLC_XY5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv"
	params:
		samplenumber = 3 #number of samples you're doing this on
	script:
		"code/utilities/germlineDEGs.R"

#upset plot of the overlap between germline DEGs in the brain and EpilC RNAseq datasets
rule Brain_EpiLC_Upset:
	input:
		"data/processed/germGENES20.csv",
		"results/DESeq2/DEGs_EpiLC_XY5cKO.csv",
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv"
	output:
		"results/figure_pieces/Upset_EpiLCBrain.pdf"
	script:
		"code/Upset_Brain_EpiLC.R"

#plot the gene ontology of EpiLC and brain DEGs
rule GO_EpiLC_vs_Brain:
	input:
		"results/DESeq2/germDEGs/germDEGs_EpiLC_XY5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv"
	output:
		"results/GO_EpiLC_vs_Brain.csv",
		"results/figure_pieces/GO_EpiLC_vs_Brain.pdf"
	script:
		"code/GO_Compare_BrainEpiLC.R"

#plot the expression of primordial germ cell markers in EpiLCs
rule EpiLC_PGC:
	input:
		"data/raw/EpiLC_gene_TPM.txt",
		"data/raw/SampleInfo_EpiLC.csv"
	output:
		"results/figure_pieces/EpiLC_PGCmarkers_XY_supplement.pdf",
		"results/figure_pieces/EpiLC_PGC_2Cell_markers_XY.pdf"
	script:
		"code/EpiLC_PGCgenes.R"

###################### Figure 4: ChIPseq KDM5C EpiLC vs PNC ############################
rule KDM5C_chip:
	input:
		"data/raw/ChIPseq_mm9_EpiLC_KDM5C_WTnoKO.bed",
		"data/raw/ChIPseq_mm9_PNC_KDM5C_WTnoKO.bed",
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_EpiLC_XY5cKO.csv"
	output:
		"data/processed/KDM5C_ChIPseq_boundpromoters_EpiLC.csv",
		"data/processed/KDM5C_ChIPseq_boundpromoters_PNC.csv",
		"results/figure_pieces/KDM5C_ChIPseq_peaklocation.pdf",
		"results/figure_pieces/KDM5C_ChIPseq_boundgermDEGs.pdf"
	script:
		"code/ChIPseq_KDM5C_EpiLC_PNC.R"

rule KDM5C_chip_GO:
	input:
		"data/processed/KDM5C_ChIPseq_boundpromoters_EpiLC.csv",
		"data/processed/KDM5C_ChIPseq_boundpromoters_PNC.csv"
	output:
		"results/ChIPseq_KDM5C_GO_EpiLC_vs_PNC.csv",
		"results/figure_pieces/KDM5C_ChIPseq_GO_all.pdf",
		"results/ChIPseq_KDM5C_GO_EpiLC_vs_PNC_unique.csv",
		"results/figure_pieces/KDM5C_ChIPseq_GO_unique.pdf",
		"results/figure_pieces/KDM5C_ChIPseq_promo_overlap.pdf"
	script:
		"code/ChIPseq_KDM5C_GO_Compare_BrainEpiLC.R"


###################### Render the manuscript ############################
rule write_paper:
	script:
		"code/utilities/paper_render.R"
