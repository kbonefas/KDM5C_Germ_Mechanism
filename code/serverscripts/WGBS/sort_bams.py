#24.08.02 - sort bam dedup alignment  

import sys
import os
import glob


#files with the deduplicated bams
INDIR = "../data/processed/BAM_dedup/"

#get the fastq files for alignment (just the first read)
BAMS = glob.glob("{INDIR}*.multiple.deduplicated.bam".format(INDIR = INDIR))

for i in BAMS:
	#clean up the read name for getting the second read
	clean = i.replace(INDIR, "")
	clean = clean.replace(".multiple.deduplicated.bam", "")

	#sort the bam file
	os.system("nice -10 samtools sort {INBAM} -o {INDIR}{clean}.sorted.dedup.bam".format(INBAM = i, INDIR = INDIR, clean = clean))
	
	#index the bam file
	os.system("nice -10 samtools index {INDIR}{clean}.sorted.dedup.bam".format(INDIR = INDIR, clean = clean))
