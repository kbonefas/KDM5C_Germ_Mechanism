#2024.08.19 Plotting average WGBS methylation at germline loci using bedgraph, divided by KDM5C bound and unbound genes
#activate the deeptools environment in conda

#convert bedgraph to bigwig using UCSC command - https://anaconda.org/bioconda/ucsc-bedgraphtobigwig
	#https://genome.ucsc.edu/goldenPath/help/bigWig.html
import os
import sys
import glob

#convert bedraphs files into bigwig files
	# get list of bedgraph files
bedgraph = glob.glob("../results/bedgraph/*.deduplicated.bedGraph.gz")
	#unzip the bedgraph file

for i in bedgraph:
	#name of file without extensions
	clean = i.replace(".bedGraph.gz", "")

	print("unzipping bedgraph for {sample}".format(sample = clean))
	# #unzip the file
	os.system("nice gzip -d {file}".format(file = i))

unzipped = glob.glob("../results/bedgraph/*.deduplicated.bedGraph")

for u in unzipped:
	clean = u.replace(".bedGraph", "")
	#make the bigwig
	print("making bigwig for {sample}".format(sample = clean))
	os.system("nice bedGraphToBigWig {inbedgraph} /nfs/value/siwase2/GENOMES/MM10/mm10.chrom.sizes {clean}.bw".format(inbedgraph = u, clean = clean))
	

