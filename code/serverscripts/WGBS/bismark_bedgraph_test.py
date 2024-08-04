#24.07.24 - bismark BAM alignment of trimmed reads. 

import sys
import os
import glob


#BAM directory
BAMDIR = "../data/processed/BAM_dedup"

#output directory
OUTDIR = "../results/bedgraph"

#get the bamfiles for bedgraph conversion
BAMS = glob.glob("{BAMDIR}/*.multiple.deduplicated.bam".format(BAMDIR = BAMDIR))

for i in BAMS:
	
	#run bismark
	os.system("nice -10 bismark_methylation_extractor --parallel 6 -o {OUTDIR} --bedGraph {bamfile}".format(OUTDIR = OUTDIR, bamfile = i))
