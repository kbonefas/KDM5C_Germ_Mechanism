#24.07.31 - bismark deduplication of BAM files

import sys
import os
import glob



#name of the samples
SAMPLES = ['KO_A_96EpiLC', 'KO_B_96EpiLC', 'WT_A_96EpiLC', 'WT_B_96EpiLC', 'KO_A_esc', 'KO_B_esc', 'WT_A_esc', 'WT_B_esc']

#output directory
OUTDIR = "../data/processed/BAM_dedup"


for i in SAMPLES:
	#get all the bam files for each sample
	BAMFILES = glob.glob("../data/processed/BAM/{sample}*.bam".format(sample = i))
	#join them together with spaces in between
	BAMFILES = " ".join(BAMFILES)

	#run bismark deduplication
	#print("{outname} - {bamfiles}".format(outname = i, bamfiles = BAMFILES))
	os.system("nice -10 deduplicate_bismark -o {outname} --output_dir {outdir} --multiple {bamfiles}".format(outname = i, outdir = OUTDIR, bamfiles = BAMFILES))
