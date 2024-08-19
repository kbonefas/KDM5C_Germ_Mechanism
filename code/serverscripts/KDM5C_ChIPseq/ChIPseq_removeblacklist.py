#remove black listed region with bedtools
#updated 23.10.31 to make it more robust
#download black listed region (BED format) at https://github.com/Boyle-Lab/Blacklist/

##Python script
import sys
import os
import glob

#where your narrowPeak files are located
PEAKpath = "../data/processed/MACS2/"
#where your blacklisted reigons bed file is
BLACKLIST = "/nfs/value/siwase2/BONEFAS/generalscripts/mm10_blacklist.bed"

samples =  glob.glob(PEAKpath + "*.narrowPeak")

for s in samples:
	#run the bedtools command
	os.system("bedtools intersect -a {sample} -b {blacklist} -v > {sample}.filtered.bed ".format(sample = s, blacklist = BLACKLIST))


