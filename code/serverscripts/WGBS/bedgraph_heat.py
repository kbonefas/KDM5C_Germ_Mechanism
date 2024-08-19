#2024.08.19 Plotting average WGBS methylation at germline loci using bedgraph, divided by KDM5C bound and unbound genes
#activate the bedtools environment in conda

import os
import sys
import glob

#bigwigs
EpiLC = ['/scratch/webserver/siwase/Testis_genes/KDM5C_H3K4me2_EpiLCs_mm10/SMCX_WT_Merged_IP_mm10.bw','/scratch/webserver/siwase/Testis_genes/KDM5C_H3K4me2_EpiLCs_mm10/SMCX_KO_Merged_IP_mm10.bw']
PNC = ['/scratch/webserver/siwase/Testis_genes/Iwase_et_al_2016_Cell_Reports_mm10/KDM5C_ChIP_WT_CN_rep1_mm10.bw','/scratch/webserver/siwase/Testis_genes/Iwase_et_al_2016_Cell_Reports_mm10/KDM5C_ChIP_KO_CN_rep1_mm10.bw']



#bedfile of the region that you want to plot
germregions = glob.glob("../data/processed/GermTSSandTES/TSS_window_3000*.bed")

#wrap all the deeptools fucntions together
#regions - gene coordinates
#bigwigs - bigwigs you're plotting
#colors - color for heatmap
def deeptoolswrap(regions, bigwigs, name):
	for i in regions:

		#name of the output - ID is the bigwig id, coord is the bed cooridnates
		coord = i.replace("../data/processed/GermTSSandTES/TSS_window_", "")
		coord = coord.replace(".bed", "")
		
		window = 3000
		out = "../results/plotProfile_TSS_KDM5C_{NAME}_{COORD}".format(COORD = coord, NAME = name)
		
		#make the compute matrix
		os.system('computeMatrix reference-point -S {BIGWIGS} -R {REGIONS} --referencePoint center -a {WINDOW} -b {WINDOW} -o {OUT}_matrix.mat.gz'.format(BIGWIGS = " ".join(bigwigs), REGIONS = i, WINDOW = window, OUT = out))
		title = "KDM5C_{COORD}".format(COORD = coord)
		
		#make the average plot --perGroup splits by bedfile so all bigwigs are plotted on same graph
		os.system('plotProfile -m {OUT}_matrix.mat.gz -out {OUT}_profile.pdf --perGroup --refPointLabel TSS --samplesLabel Kdm5c_WT Kdm5c_KO --plotTitle {TITLE} --plotWidth 10 --plotHeight 10 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#make heatmap plot
		os.system('plotHeatmap -m {OUT}_matrix.mat.gz -out {OUT}_heatmap.pdf --perGroup --yMax 40 --zMax 60 --colorMap RdPu --refPointLabel TSS --missingDataColor white --samplesLabel Kdm5c_WT Kdm5c_KO --plotTitle {TITLE} --heatmapHeight 14 --heatmapWidth 5 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#os.system("plotHeatmap -m {OUT}_matrix.mat.gz --colorList 'white, #63344c' 'white, #63344c' 'white, #275c62' 'white, #275c62' --missingDataColor white --heatmapHeight 14 --heatmapWidth 5 -out {OUT}_heatmap.pdf".format(MATRIX = MATRIX, NAME = NAME))


deeptoolswrap(germregions, EpiLC, 'EpiLC')
deeptoolswrap(germregions, PNC, 'PNC')