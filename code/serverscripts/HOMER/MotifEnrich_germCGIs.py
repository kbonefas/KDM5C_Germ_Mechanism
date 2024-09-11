# 24.07.19 HOMER to find motifs enriched in KDM5C bound vs unbound germline DEGs to see what is promoting their dysregulation independent of KDM5C

import sys
import os
import glob

#where HOMER is located
HOMER = "../../HOMER/"

#output directory
OUT = "../results"

#window of bp that you want to search aroudn the TSS
WINDOW = "500"

#input files of gene list
#input file needs to be a .txt file with the gene IDs in the first column, and no quotes around them.
GENES = glob.glob("../data/raw/*CGI_HOMER.txt")

for i in GENES:
	#clean the file name 
	clean = i.replace("../data/raw/", "")
	clean = clean.replace("_HOMER.txt", "")

	#make a directory for the outputs for that sample because for some reason that's what homer does
	OUTDIR = "{OUT}/{CLEAN}".format(OUT = OUT, CLEAN = clean)
	os.system("mkdir {OUTDIR}".format(OUTDIR = OUTDIR))

	#run the HOMER command for that gene list
	os.system("nohup perl {HOMER}bin/findMotifs.pl {GENES} mouse {OUTDIR} -start -{WINDOW} -end {WINDOW} &".format(HOMER = HOMER, OUT = OUT, GENES = i, OUTDIR = OUTDIR, WINDOW = WINDOW))
