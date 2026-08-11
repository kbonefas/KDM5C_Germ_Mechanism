#peakcalling with MACS2 for KDM5C ESC peaks 
#2026.08.11

#before running the command check the follwing. 
#is macs2 installed and your conda environment activated?
#Is your input data deep enough? Combine all the input bam files. 50M reads are needed at least. Deeper the input, more precise the peak calling. 

import sys
import os
import glob

#where your bam files are located
BAMpath = "/nfs/turbo/umms-siwase/SHIGEKI2/kdm5c_germline_genes/data/kdm5c/bam/"


#macs2 flag options
#format is BAM (if using paired end reads use BAMPE), genome is mouse (mm), q cutoff, verbose level 3 to get debug info
MACS2_FLAGS = "-f BAM -g mm -q 0.05 --verbose 3"
INPUT = BAMpath + "sort_Agarwal_Input_merged.bam"

samples = glob.glob(BAMpath + "sort_Agarwal_KDM5C_*.bam")

for s in samples:
    clean = s.replace(BAMpath, "")
    clean = clean.replace(".bam", "")
    clean = clean.replace("sort_", "")

    os.system("macs2 callpeak {FLAGS} -t {SAMPLE} -c {INPUT} -n {clean} --outdir /nfs/turbo/umms-siwase/SHIGEKI2/kdm5c_germline_genes/data/kdm5c/macs2".format(FLAGS = MACS2_FLAGS, SAMPLE = s, INPUT = INPUT, clean = clean))


