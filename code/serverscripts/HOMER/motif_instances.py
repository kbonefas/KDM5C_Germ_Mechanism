#2024.07.19 - find specific instances of motifs using HOMER
#activate blueberry conda env to get python 3
#have to run without nohup to avoid the running output being put on the top of the file

#24.07.23 - using my custom e2f and ebox motifs

import sys
import os
import glob

#where HOMER is located
HOMER = "../../HOMER/"

#output directory
OUT = "../results"

#make a dictionary of the custom motifs
motifs = {"E2F": "TCCCGC", "Ebox": "CACGTG", "xbox" : "usingHOMERs"}


#make a directory for the outputs for that sample because for some reason that's what homer does
outdir = "{OUT}/MotifsCount".format(OUT = OUT)
#os.system("mkdir {OUTDIR}".format(OUTDIR = outdir))



#function for searching for motifs in list of genes
#MOTIF - motif you're searching for
#GENES - genes you're searching (txt file)
#OUTDIR - output directory
#CLEAN - cleaned file name
def motifSearch(MOTIF, GENES):

	#set a header so x-box motifs can use the homer provided one
	if MOTIF == "xbox":
		header = "{HOMER}motifs/".format(HOMER = HOMER)
	else:
		header = ""
		
	os.system("perl {HOMER}bin/findMotifs.pl {GENES} mouse {OUTDIR}/ -find {header}{MOTIF}.motif > {OUTDIR}/{MOTIF}_instances_findMotifs_{CLEAN}.txt".format(HOMER = HOMER, MOTIF = MOTIF, header = header, GENES = GENES, OUTDIR = outdir, CLEAN = clean))

	#os.system("perl {HOMER}bin/annotatePeaks.pl tss mm10 -list {GENES} -m {header}{MOTIF}.motif -size -500,500 > {OUTDIR}/{MOTIF}_instances_annoPeaks_300_300_{CLEAN}.txt".format(HOMER = HOMER, MOTIF = MOTIF, GENES = GENES, header = header, OUTDIR = outdir, CLEAN = clean))
	
	#histogram of motif instances
	#os.system("perl {HOMER}bin/annotatePeaks.pl tss mm10 -list {GENES} -m {MOTIF}.motif -size -500,250 -hist 10 > {OUTDIR}/{MOTIF}_hist_{CLEAN}.txt".format(HOMER = HOMER, MOTIF = MOTIF, GENES = GENES, OUTDIR = outdir, CLEAN = clean))




#get input files of gene list
#input file needs to be a .txt file with the gene IDs in the first column, and no quotes around them.
GENES = glob.glob("../data/raw/*.txt")


#for every motif in the dictionary, make the motif file then run matches in each gene list
for x, y in motifs.items():
	if x == "E2F" or x == "Ebox" :
		#make the custom motif file
		os.system("{HOMER}bin/seq2profile.pl {SEQ} 0 {NAME} > {NAME}.motif".format(HOMER = HOMER, SEQ = y, NAME = x))
		
		#search for motifs in gene list		
		for i in GENES:
			#clean the file name 
			clean = i.replace("../data/raw/", "")
			clean = clean.replace("_HOMER.txt", "")
			
			motifSearch(x, i)	


	else:
	#search for motifs in gene list
		for i in GENES:
			#clean the file name 
			clean = i.replace("../data/raw/", "")
			clean = clean.replace("_HOMER.txt", "")
			
			motifSearch(x, i)	




