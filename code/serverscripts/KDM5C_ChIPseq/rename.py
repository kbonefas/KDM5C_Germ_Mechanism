## script to rename the BAM names from core to the actual sample names ##
import os
import glob

# get list of samples
#folder with files
filefolder = "../data/processed/BAM/"
samples =  glob.glob(filefolder + "*.nodup.bam")
	
	#for each sample, rename the sample based on the matching key
for s in samples:
	#print(s)
	#clean up the sample names by removing extraneous info
	idname1 = s.replace(filefolder, "")
	idname = idname1.replace(".mm10.nodup.bam", "")

	#open the key file
	key = open("../data/raw/KDM5C_EpiLC_ChIP_sampleinfo.txt", "r")
	for line in key:
		#split the line based on the tab delimiter
		match = line.split()

		#if the sample ID matches the one in the file
		if idname[0:len(idname)] == match[0]:
			#write what the name is being changed to
			print(match[0] + " = " + match[1])
			#rename the files
			os.system("mv {name} {filefolder}{newname}.mm10.bam".format(name = s, filefolder = filefolder, newname = match[1]))
			#rename the bam.bai files
			os.system("mv {filefolder}{idname1}.bai {filefolder}{newname}.mm10.bam.bai".format(idname1= idname1, filefolder = filefolder, newname = match[1]))
			#close the key file
			key.close() 
			break;








