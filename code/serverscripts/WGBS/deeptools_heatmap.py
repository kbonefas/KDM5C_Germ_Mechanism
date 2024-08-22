#2024.08.19 Plotting average WGBS methylation at germline loci using bedgraph, divided by KDM5C bound and unbound genes
#activate the deeptools environment in conda

#convert bedgraph to bigwig - https://anaconda.org/bioconda/ucsc-bedgraphtobigwig

import os
import sys
import glob

#bigwigs of EpiLCs
#combine the bigwigs together to make an average
bigwigs = glob.glob("../results/bedgraph/*EpiLC*.bw")

genotypes = ['WT', 'KO']

for g in genotypes:
	bw = [k for k in bigwigs if g in k]

	#make an average bigwig
	os.system("bigwigAverage -b bw[0] bw[1] -o ../results/bedgraph/avg_exEpiLC_{genotype}.bw".format(genotype = g))

#one graph
	#WT and 5CKO bigwigs
	#regions - CGI that are KDM5C bound and unbound?


#bedfile of the region that you want to plot
avgbigwigs = glob.glob("../results/bedgraph/avg_exEpiLC*.bw")

#bed file coordinates
germ_regions = ["TSS", "CGI"]

#wrap all the deeptools fucntions together
#regions - gene coordinates
#bigwigs - bigwigs you're plotting
#colors - color for heatmap
def deeptoolswrap(regions, bigwigs):
	for i in regions:

		#get the bedfiles of that region
		beds = glob.glob("../data/processed/*{reg}*.bed".format(reg = i))
		beds = " ".join(beds)

		#name of the output 
		out = "../results/deeptools/plotProfile_WGBS_EpiLC_WTvsKO_{reg}".format(reg = i)
		
		#make the compute matrix
		os.system('computeMatrix scale-regions -S {BIGWIGS} -R {REGIONS} -o {OUT}_matrix.mat.gz'.format(BIGWIGS = " ".join(bigwigs), REGIONS = beds, OUT = out))
		title = "Germline_{reg}".format(reg = i)
		
		#make the average plot --perGroup splits by bedfile so all bigwigs are plotted on same graph
		os.system('plotProfile -m {OUT}_matrix.mat.gz -out {OUT}_profile.pdf --perGroup --startLabel start --endLabel end --samplesLabel Kdm5c_KO Kdm5c_WT --plotTitle {TITLE} --plotWidth 10 --plotHeight 10 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#make heatmap plot
		#--yMax 40 --zMax 60
		os.system('plotHeatmap -m {OUT}_matrix.mat.gz -out {OUT}_heatmap.pdf --perGroup --colorMap RdPu --startLabel start --endLabel end --missingDataColor white --samplesLabel Kdm5c_KO Kdm5c_WT --plotTitle {TITLE} --heatmapHeight 14 --heatmapWidth 5 --plotFileFormat pdf'.format(OUT = out, TITLE = title))
		
		#os.system("plotHeatmap -m {OUT}_matrix.mat.gz --colorList 'white, #63344c' 'white, #63344c' 'white, #275c62' 'white, #275c62' --missingDataColor white --heatmapHeight 14 --heatmapWidth 5 -out {OUT}_heatmap.pdf".format(MATRIX = MATRIX, NAME = NAME))


deeptoolswrap(germ_regions, avgbigwigs)
