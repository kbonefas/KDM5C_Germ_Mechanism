#fastqc analysis of reads. Make sure you have fastqc installed

import sys
import os
import glob


#fastq file input location
FASTQ = "../data/raw"

#output directory
OUTDIR = "../data/processed/FASTQC"


os.system("fastqc -o {outdir} -t 2 {fastq}/*.fq.gz".format(outdir = OUTDIR, fastq = FASTQ))