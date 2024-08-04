#24.07.24 - bismark BAM alignment of trimmed reads. 

import sys
import os
import glob


#fastq file input location (trimmed)
FASTQ_DIR = "../data/processed/trimmed_reads/"

#output directory
OUTDIR = "../data/processed/BAM"

#get the fastq files for alignment (just the first read)
firstREAD = glob.glob("{fastq_dir}*val_1.fq.gz".format(fastq_dir = FASTQ_DIR))

for i in firstREAD:
	#clean up the read name for getting the second read
	clean = i.replace(FASTQ_DIR, "")
	clean = clean.replace("1_val_1.fq.gz", "")

	#run bismark
	os.system("nice -10 bismark --non_directional --parallel 6 -o {outdir} /nfs/value/siwase2/GENOMES/MM10/bismark_genome/ -1 {READ_1} -2 {fastq_dir}{READ_2}2_val_2.fq.gz".format(outdir = OUTDIR, READ_1 = i, fastq_dir = FASTQ_DIR, READ_2 = clean))
