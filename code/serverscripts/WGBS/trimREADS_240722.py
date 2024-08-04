#fastqc analysis of reads. Make sure you have fastqc installed
#activate bismark conda environment

import sys
import os
import glob


#fastq file input location
FASTQ = "../data/raw"

#output directory
OUTDIR = "../data/processed/trimmed_reads"


os.system("trim_galore -o {outdir} --paired {fastq}/*.fq.gz".format(outdir = OUTDIR, fastq = FASTQ))
