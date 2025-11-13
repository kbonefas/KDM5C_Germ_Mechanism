##### 2025.11.12
### plot heatmaps of previously published KDM5C and PRC1.6 ChIP-seq in ESCs
import sys
import os
import glob



KDM5C1 = ['Bigwig_Analysis/Agarwal_KDM5C_rep1.bw','Bigwig_Analysis/Agarwal_Input_rep1.bw']
KDM5C2 = ['Bigwig_Analysis/Agarwal_KDM5C_rep2.bw','Bigwig_Analysis/Agarwal_Input_rep2.bw']
MGA = ['Bigwig_Analysis/Stielow_MGA.bw','Bigwig_Analysis/Stielow_cntrlIgG.bw']
PCGF6 = ['Bigwig_Analysis/Stielow_PCGF6.bw','Bigwig_Analysis/Stielow_cntrlIgG.bw']



#bedfile of the region that you want to plot
germregions = ["TSS_window_3000bp_all_germ.bed"]

#wrap all the deeptools fucntions together
#regions - gene coordinates
#bigwigs - bigwigs you're plotting
#colors - color for heatmap
def deeptoolswrap(regions, bigwigs, name):
	for i in regions:

		#name of the output - ID is the bigwig id, coord is the bed cooridnates
		coord = "all_germ"
		
		window = 3000
		out = "plotProfile_TSS_KDM5C_{NAME}_{COORD}".format(COORD = coord, NAME = name)
		
		#make the compute matrix
		os.system('computeMatrix reference-point -S {BIGWIGS} -R {REGIONS} --referencePoint center -a {WINDOW} -b {WINDOW} -o {OUT}_matrix.mat.gz'.format(BIGWIGS = " ".join(bigwigs), REGIONS = i, WINDOW = window, OUT = out))
		title = "{NAME}_{COORD}".format(NAME = name, COORD = coord)
		
		#make the average plot --perGroup splits by bedfile so all bigwigs are plotted on same graph
		os.system('plotProfile -m {OUT}_matrix.mat.gz -out {OUT}_profile.pdf --perGroup --refPointLabel TSS --samplesLabel Kdm5c_WT Kdm5c_KO --plotTitle {TITLE} --plotWidth 10 --plotHeight 10 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#make heatmap plot
		os.system('plotHeatmap -m {OUT}_matrix.mat.gz -out {OUT}_heatmap.pdf --perGroup --yMax 40 --zMax 60 --colorMap RdPu --refPointLabel TSS --missingDataColor white --samplesLabel Kdm5c_WT Kdm5c_KO --plotTitle {TITLE} --heatmapHeight 14 --heatmapWidth 5 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#os.system("plotHeatmap -m {OUT}_matrix.mat.gz --colorList 'white, #63344c' 'white, #63344c' 'white, #275c62' 'white, #275c62' --missingDataColor white --heatmapHeight 14 --heatmapWidth 5 -out {OUT}_heatmap.pdf".format(MATRIX = MATRIX, NAME = NAME))


deeptoolswrap(germregions, KDM5C1, 'KDM5C1')
deeptoolswrap(germregions, KDM5C2, 'KDM5C2')
deeptoolswrap(germregions, MGA, 'MGA')
deeptoolswrap(germregions, PCGF6, 'PCGF6')