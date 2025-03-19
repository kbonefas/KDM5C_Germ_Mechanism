
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



###################### Germline-enriched Genes ############################
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
		"data/raw/20150917_Soh_embryonicmouseWWvsequencing.xlsx", #FPKM of embryonic WWv (germline-depeleted) and WT germline
			#data from https://pubmed.ncbi.nlm.nih.gov/26378784/
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


###################### EpiLC vs Brain germline genes ############################

#DESeq2 on male and female Kdm5c mutant EpiLCs
rule EpiDESeq2_XXvsXY:
	input:
		cts = "data/raw/EpiLC_gene_expected_count.txt"
	output:
		"data/raw/SampleInfo_EpiLC_XXvsXY.csv",
		"results/figure_pieces/PCA_EpiLC.pdf",
		"data/processed/restable_EpiLC_XXvsXY_XY5cKO.csv",
		"data/processed/restable_EpiLC_XXvsXY_XX5cHET.csv",
		"data/processed/restable_EpiLC_XXvsXY_XX5cKO.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cHET.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cKO.csv"
	params:
		alpha = PADJ
	script:
		"code/DESeq2_EpiLC_XXvsXY.R"

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
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
		"data/processed/restable_EpiLC_XXvsXY_XY5cKO.csv"

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
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
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
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
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
		"results/figure_pieces/EpiLC_PGC_2Cell_markers_XY.pdf",
		"results/figure_pieces/EpiLC_piRNAgenes_TPM_XY.pdf"		
	script:
		"code/EpiLC_PGCgenes.R"


################ male vs female EpiLCs ################

#plot the overlap between XX and XY EpiLC 5cKO/HET DEGs and which chromosome
rule EpiLC_XXvsXY:
	input:
		"data/processed/germGENES20.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cHET.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cKO.csv"
	output:
		"results/EpiLC_XXvsXY_germDEGs.csv",
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_euler.pdf",
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_shared_chromo.pdf", 
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_XXonly_chromo.pdf", 
		"results/EpiLC_XXvsXY_XXDEGs_chromosome.csv", #How many XX-specific DEGs are on each chromosome 
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_histogram.pdf",
		"results/EpiLC_XXvsXY_germlineDEGs_GO.csv",
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_GO.pdf",
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_allXXDEGs_chromo.pdf",
		"results/figure_pieces/EpiLC_XXvsXY_germlineDEGs_chromo_proportion.pdf"
	script:
		"code/EpiLC_XXvsXY.R"

rule EpiLC_XXvsXY_heat:
	input:
		"data/processed/EpiLC_XXvsXY_germDEGs.csv",
		"data/processed/restable_EpiLC_XXvsXY_XY5cKO.csv",
		"data/processed/restable_EpiLC_XXvsXY_XX5cHET.csv",
		"data/processed/restable_EpiLC_XXvsXY_XX5cKO.csv"
	output:
		"results/figure_pieces/EpiLC_XXvsXY_l2fc_heat.pdf",
		"data/processed/EpiLC_XXvsXY_heat_l2fc_clusters.csv",
		"results/figure_pieces/EpiLC_XXvsXY_l2fc_heat_clustnumb.pdf"
	script:
		"code/EpiLC_XXvsXY_heat.R"

################ Egg vs sperm genes in EpiLCs ################

#get which germline genes are egg/sperm-biased and plot in EpiLCs
rule eggvssperm:
	input:
		"data/processed/germGENES20.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cHET.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XX5cKO.csv"
	output:
		"results/figure_pieces/GermGenes_eggvssperm_sankey.pdf",
		"results/figure_pieces/GermGenes_eggvssperm_XXvsXYEpiLC.pdf"
	script:
		"code/GermGenes_eggvssperm.R"


###################### Figure 5: ChIPseq KDM5C EpiLC vs PNC ############################
rule KDM5C_chip:
	input:
		"data/raw/ChIPseq_mm10_EpiLC_WTnoKO_consensus_peaks.bed",
		"data/raw/ChIPseq_mm10_PNC_WTnoKO_consensus_peaks.bed",
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_EpiLC_XY5cKO.csv",
		"data/processed/germGENES20.csv"
	output:
		"data/processed/KDM5C_ChIPseq_boundpromoters_EpiLC.csv",
		"data/processed/KDM5C_ChIPseq_boundpromoters_PNC.csv",
		"results/figure_pieces/KDM5C_ChIPseq_peaklocation.pdf",
		"results/figure_pieces/KDM5C_ChIPseq_boundgermDEGs.pdf",
		"results/KDM5C_binding_germDEGs_EpiLC.csv",
		"data/processed/KDM5C_bound_germDEGS_HOMER.txt",
		"data/processed/KDM5C_unbound_germDEGs_HOMER.txt",
		"results/KDM5C_binding_allgerm_EpiLC.csv",
		"data/processed/KDM5C_bound_allgerm_HOMER.txt",
		"data/processed/KDM5C_unbound_allgerm_HOMER.txt",
		"results/figure_pieces/KDM5C_ChIPseq_germ_euler.pdf"
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

rule KDM5C_chip_HOMER:
	input:
		#all of the kdm5c bound and unbound germ degs
		"results/KDM5C_binding_allgerm_EpiLC.csv",

		#germline genes with Ebox and E2f motifs near TSS
		"data/raw/HOMER/E2F_instances_findMotifs_KDM5C_bound_allgerm.txt",
		"data/raw/HOMER/Ebox_instances_findMotifs_KDM5C_bound_allgerm.txt",
		"data/raw/HOMER/E2F_instances_findMotifs_KDM5C_unbound_allgerm.txt",
		"data/raw/HOMER/Ebox_instances_findMotifs_KDM5C_unbound_allgerm.txt",
		"data/raw/HOMER/xbox_instances_findMotifs_KDM5C_bound_allgerm.txt",
		"data/raw/HOMER/xbox_instances_findMotifs_KDM5C_unbound_allgerm.txt"

	output:
		"results/figure_pieces/KDM5C_ChIPseq_HOMER_e2febox_perc_bar.pdf",
		"results/figure_pieces/KDM5C_ChIPseq_HOMER_xbox_perc_bar.pdf"
	script:
		"code/ChIPseq_KDM5C_HOMER.R"

rule KDM5C_ChIP_Rfx2:
	input:
		"data/raw/EpiLC_gene_TPM.txt",
		"data/raw/SampleInfo_EpiLC.csv"
	output:
		"results/figure_pieces/EpiLC_Rfx2_mRNA.pdf"
	script:
		"code/EpiLC_Rfx2.R"

rule KDM5C_chip_Stra8:
	input:
		"results/KDM5C_binding_allgerm_EpiLC.csv",
		"data/raw/Kojima_eLife_Stra8_bound_genes.csv", #from https://elifesciences.org/articles/43738 
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_EpiLC_XY5cKO.csv"

	output:
		"results/figure_pieces/KDM5C_ChIPseq_Stra8_targets.pdf",
		"results/figure_pieces/KDM5C_ChIPseq_Stra8_CpG_island.pdf"
	script:
		"code/ChIPseq_KDM5C_Stra8.R"

###################### Figure 5: Germline DEGs and RA signaling ############################
germgenes = "data/processed/germGENES20.csv"
#use DESeq2 to test which genes are significantly changed between WT and 5CKO for each condition 
rule ESCEpiLC_RA_WTKO_DESeq2:
	input:
		cts = "data/raw/230919_ESCEpiLC_RA_gene_expected_count.annot.txt",
		sampleinfo = "data/raw/SampleInfo_ESCEpiLC_RA.csv",
		germ = germgenes
	output:
		"results/figure_pieces/PCA_ESCEpiLC_RA_5CKOvWT.pdf",
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_0.csv",
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_48RA.csv",		
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_48NO.csv",
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_96RA.csv",		
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_96NO.csv",
		"results/DESeq2/germDEGs/germDEGs_ESCEpiLC_RA_5CKOvWT_0.csv",
		"results/DESeq2/germDEGs/germDEGs_ESCEpiLC_RA_5CKOvWT_48RA.csv",
		"results/DESeq2/germDEGs/germDEGs_ESCEpiLC_RA_5CKOvWT_48NO.csv",
		"results/DESeq2/germDEGs/germDEGs_ESCEpiLC_RA_5CKOvWT_96RA.csv",
		"results/DESeq2/germDEGs/germDEGs_ESCEpiLC_RA_5CKOvWT_96NO.csv"
	params:
		alpha = PADJ
	script:
		"code/DESeq2_ESC_EpiLC_RA_WTvKO.R"

rule ESC_EpiLC_markers:
	input:
		"data/raw/230919_ESCEpiLC_RA_gene_TPM.annot.txt", #annotated TPM for DESeq2
		"data/raw/SampleInfo_ESCEpiLC_RA.csv"
	output:
		"results/figure_pieces/ESC_EpiLC_marker_heatmap_naive.pdf",
		"results/figure_pieces/ESC_EpiLC_marker_heatmap_primed.pdf"
	script:
		"code/ESC_EpiLC_markers.R"

#graphing of germline genes in ESC and EpiLCs
rule ESC_EpiLC_5CKOcluster:
	input:
		germgenes,
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_0.csv",	
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_48NO.csv",
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_48RA.csv",		
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_96NO.csv",
		"data/processed/restable_ESCEpiLC_RA_5CKOvWT_96RA.csv",	
		"data/raw/230919_ESCEpiLC_RA_gene_TPM.annot.txt",
		"data/raw/SampleInfo_ESCEpiLC_RA.csv",
		"results/DESeq2/germDEGs/germDEGs_hip5cKO.csv",
		"results/DESeq2/germDEGs/germDEGs_amy5cKO.csv"
	output:
		"results/figure_pieces/ESCEpiLC_RA_allgermDEGs_hist.pdf",
		"results/figure_pieces/ESCEpiLC_RA_5CKO_germ_heat_TPM.pdf",
		"results/ESCEpiLC_RA_5CKO_germ_TPM_clusters.csv",
		"results/figure_pieces/ESCEpiLC_RA_5CKO_germ_heat_TPM_clustnumb.pdf",
		"results/ESCEpiLC_RA_5CKO_germDEGs_96RAonly.csv",
		"data/processed/ESCEpiLC_RA_germ_TPM_5CKOclusters_hip5cKODEGs.csv",
		"data/processed/ESCEpiLC_RA_germ_TPM_5CKOclusters_amy5cKODEGs.csv",
		"data/processed/ESCEpiLC_RA_5CKOclusters_1.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_2.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_3.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_4.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_5.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_6.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_7.bed",
		"data/processed/ESCEpiLC_RA_5CKOclusters_8.bed"
	params:
		alpha = PADJ
	script:
		"code/ESC_EpiLC_5CKOcluster.R"

#gene ontology of 5CKO clusters
rule ESC_EpiLC_5CKOcluster_GO:
	input:
		"results/ESCEpiLC_RA_5CKO_germ_TPM_clusters.csv"
	output:
		"results/ESCEpiLC_RA_germ_TPM_5CKOclusters_GO.csv",
		"results/figure_pieces/ESCEpiLC_RA_germ_TPM_5CKOclusters_GO.pdf"
	script:
		"code/ESC_EpiLC_5CKOcluster_GO.R"


#TPM plot for RA+/- WT and 5CKO of whichever gene you input
rule ESCEpiLC_RA_TPM:
	input:
		"data/raw/230919_ESCEpiLC_RA_gene_TPM.annot.txt",
		"data/raw/SampleInfo_ESCEpiLC_RA.csv"
	output:
		"results/figure_pieces/ESCEpiLC_RA_germ_5CKO_RAgenesTPM.pdf",
		"results/figure_pieces/ESCEpiLC_RA_germ_5CKOonly_RAgenesTPM.pdf",
		"results/figure_pieces/ESCEpiLC_RA_TPM_genesofinterest.pdf"
	script:
		"code/ESC_EpiLC_RA_TPM.R"


rule RA_ICC:
	input:
		"data/raw/Cell_Counts/Dazl_100nMRA/241009_Dazl_100nMRA.xlsx"
	output:
		"results/figure_pieces/RA_ICC_DAZL_100nMRA.pdf"
	script:
		"code/RA_ICCplots.R"


############# Tissue-enriched DEGs #############

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

rule tissue_genes_dot:
	input:
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv",
	output:
		"results/figure_pieces/TissueSpecific_dot.pdf"
	script:
		"code/Tissue_genes_dot.R"



###################### KDM5C substrates, CGIs, and DNAme ############################
rule GeneTSS:
	input:
		"data/processed/germGENES20.csv",
		"results/DESeq2/DEGs_EpiLC_XXvsXY_XY5cKO.csv",
		"results/DESeq2/DEGs_amy5cKO.csv",
		"results/DESeq2/DEGs_hip5cKO.csv"
	output:
		"data/processed/TSS_to_TES_all_germ.bed",
		"data/processed/TSS_window_3000bp_all_germ.bed",		
		"data/processed/TSS_to_TES_testgerm.bed",
		"data/processed/TSS_window_3000bp_testgerm.bed",
		"data/processed/TSS_to_TES_all_germ_DEGs.bed",
		"data/processed/TSS_window_3000bp_all_germ_DEGs.bed",
		"data/processed/TSS_to_TES_nongerm_upDEGs.bed",
		"data/processed/TSS_window_3000bp_nongerm_upDEGs.bed"
	script:
		"code/GeneTSS.R"

rule KDM5C_ESCtoEpiLC:
	input:
		"data/raw/220303_Kdm5c_ESCEpiLC_western_quant.csv", #quantification spreadsheet
		"data/raw/220406_ESCtoEpiLC_Kdm5c_qPCR.csv"
	output:
		"results/figure_pieces/KDM5C_ESCtoEpilC_western.pdf",
		"results/figure_pieces/KDM5C_ESCtoEpilC_qPCR.pdf"
	params:
		expprotein = "KDM5C", #experimental protein
		housekeep = "DAXX" #housekeeping gene
	script:
		"code/KDM5C_ESCtoEpiLC.R"

############## KDM5C WGBS ######################

#make the bed file for calculating methylation changes
rule WGBS_germ_regions:
	input:
		"data/processed/germGENES20.csv",
		"data/raw/EpiLC_gene_expected_count.txt"
	output:
		"data/processed/TSS_window_500bp_all_germ.bed",
		"data/processed/TSS_window_500bp_all_genes.bed"
	script:
		"code/WGBS_germ_regions.R"

### CGI germline genes ###
#Find how many germline genes have CGIs at promoter
rule Germ_CGI:
	input:
		"data/raw/CGI_UCSC.bed",
		"results/KDM5C_binding_allgerm_EpiLC.csv"
	output:
		"data/processed/CGI_UCSC_all_germ.bed",
		"results/KDM5C_binding_allgerm_CGI.csv",
		"results/figure_pieces/WGBS_CGI_bar.pdf",
		"results/germline_CGI_GO.csv",
		"results/figure_pieces/WGBS_CGI_GO.pdf",
		"data/processed/germ_CGI_HOMER.txt",
		"data/processed/germ_noCGI_HOMER.txt",
	script:
		"code/WGBS_germ_CGI.R"

#stage of germ cell development CGI vs non-CGI genes are expressed 
rule Germ_CGI_stage:
	input:
		"results/KDM5C_binding_allgerm_CGI.csv",
		"data/raw/Green_2018_logAvgNormalizedExpression_GermCell.csv"
			#from green 2018 (https://pubmed.ncbi.nlm.nih.gov/30146481/)
				# GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp
	output:
		"results/figure_pieces/CGIvsNON_germcellstages.pdf",
		"results/figure_pieces/CGIvsNON_germcellstages_heat_CGI.pdf",
		"results/figure_pieces/CGIvsNON_germcellstages_heat_noCGI.pdf"
	script:
		"code/Germ_CGI_stage.R"


#Get the genes for each location and make volcano plots
rule WGBS_volcano:
	input:
		"results/KDM5C_binding_allgerm_CGI.csv",
		"data/processed/WGBS_restab_regionCounts_allgermTSS500_WT_min3_ESCvsEpiLC.csv",
		"data/processed/WGBS_restab_regionCounts_allgermCGI_WT_min3_ESCvsEpiLC.csv",
		"data/processed/WGBS_restab_regionCounts_allgermTSS500_min3_WTvsKO_EpiLC.csv",
		"data/processed/WGBS_restab_regionCounts_allgermCGI_min3_WTvsKO_EpiLC.csv",
		"data/processed/WGBS_restab_regionCounts_allgenesTSS500_min3_WTvsKO_EpiLC.csv",
		"data/raw/EpiLC_gene_expected_count.txt"

	output:
		"results/WGBS/WGBS_restab_annotated_allgermTSS_500bp_WT_ESCvsEpiLC.csv",
		"results/WGBS/WGBS_restab_annotated_allgermCGI_WT_ESCvsEpiLC.csv",
		"results/figure_pieces/WGBS_volcano_ESCvsEpiLC_allgermTSS500.pdf",
		"results/figure_pieces/WGBS_volcano_ESCvsEpiLC_allgermCGI.pdf",
		"results/WGBS/WGBS_restab_annotated_allgermTSS_500bp_WTvsKO.csv",
		"results/WGBS/WGBS_restab_annotated_allgermCGI_WTvsKO.csv",
		"results/figure_pieces/WGBS_volcano_EpiLC_WTvsKO_allgermTSS500.pdf",
		"results/figure_pieces/WGBS_volcano_EpiLC_WTvsKO_allgermCGI.pdf",		
		"results/WGBS/WGBS_restab_annotated_allgenes_WTvsKO.csv",
		"results/figure_pieces/WGBS_volcano_EpiLC_WTvsKO_allgenesTSS500_germ.pdf",
		"results/figure_pieces/WGBS_GO_EpiLC_WTvsKO_allgenesTSS500_hypo.pdf",
		"results/GO_WGBS_exEpiLC_WTvsKO_allgenes_hypo.csv"
	script:
		"code/WGBS_volcano.R"

#histogram of percent methylation at germline promoters
rule WGBS_hist:
	input:
		"results/KDM5C_binding_allgerm_CGI.csv",
		"data/processed/WGBS_percmeth_allgermTSS500_min3.csv",
		"results/WGBS/WGBS_restab_annotated_allgermTSS_500bp_WT_ESCvsEpiLC.csv"
	output:
		"results/WGBS/WGBS_percmeth_annotated_allgermTSS_500bp.csv",
		"results/figure_pieces/WGBS_hist_percMeth_5CKO.pdf"
	script:
		"code/WGBS_germ_hist.R"


###################### Render the manuscript ############################
rule write_paper:
	script:
		"code/utilities/paper_render.R"
