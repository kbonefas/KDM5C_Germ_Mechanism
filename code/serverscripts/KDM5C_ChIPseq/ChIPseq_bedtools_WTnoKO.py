#24.05.14 - Remove the KO peaks from the WT samples using bedtools

import sys
import os
import glob

#deactivate conda and activate bluberry environment

outfolder = "../data/processed/"



#running bedtools intersect to remove overlapping KO peaks from WT peaks
def bedtools_command(samples, ID):
	WTbed = [item for item in samples if 'wt' in item.lower()][0]
	KObed = [item for item in samples if 'ko' in item.lower()][0]


	bedtools_string = "bedtools intersect -v -a {WTbed} -b {KObed} > {outfolder}{ID}_WTnoKO_consensus_peaks.bed".format(outfolder = outfolder, WTbed = WTbed, KObed = KObed, ID = ID)
	os.system(bedtools_string)

#for PNCs
PNCsamps = glob.glob("../data/processed/consensus*PNC*.bed")
EpiLCsamps = glob.glob("../data/processed/MACS2/epilc*filtered.bed")


bedtools_command(PNCsamps, "PNC")
bedtools_command(EpiLCsamps, "EpiLC")
