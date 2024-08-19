#peakcalling with MACS2
#script updated 23.10.31 by Katie Bonefas to make it more robust
#24.05.10 Update for specific case

#before running the command check the follwing. 
#is macs2 installed and your conda environment activated?
#Is your input data deep enough? Combine all the input bam files. 50M reads are needed at least. Deeper the input, more precise the peak calling. 

import sys
import os
import glob

#where your bam files are located
BAMpath = "../data/processed/BAM/"


#macs2 flag options
#format is BAM (if using paired end reads use BAMPE), genome is mouse (mm), q cutoff, verbose level 3 to get debug info
MACS2_FLAGS = "-f BAM -g mm -q 0.05 --verbose 3"
INPUT = BAMpath + "allInput.merged.bam"


#make a macs2 function for these specific files
#sample id is if it's wt or ko
def MACS2(sampleid):

	samples = glob.glob(BAMpath + "*{sampleid}*.bam".format(sampleid = sampleid))
	
	if sampleid == 'wt':
		samples = [item for item in samples if 'merged' in item]
	
		print(INPUT)
	
		for s in samples:
			#clean up the sample names
			clean = s.replace(BAMpath, "")
			clean = clean.replace(".mm10.bam", "")
	
			#skip any input samples
			clean = clean.lower()
			if 'input' in clean:
				continue
	
			#if it's not an input sample
			else:
				print("running macs2 for " + clean)
				#run the macs2 command
				#os.system("macs2 callpeak {FLAGS} -t {SAMPLE} -c {INPUT} -n {clean} --outdir ../data/processed/MACS2".format(FLAGS = MACS2_FLAGS, SAMPLE = s, INPUT = INPUT, clean = clean))
				#os.system('sleep 20')



	else:
		for s in samples:
			#clean up the sample names
			clean = s.replace(BAMpath, "")
			clean = clean.replace(".mm10.bam", "")
	
			#skip any input samples
			clean = clean.lower()
			print(clean)
			if 'input' in clean:
				continue
	
			#if it's not an input sample
			else:
				print("running macs2 for " + clean)
				#run the macs2 command
				os.system("macs2 callpeak {FLAGS} -t {SAMPLE} -c {INPUT} -n {clean} --outdir ../data/processed/MACS2".format(FLAGS = MACS2_FLAGS, SAMPLE = s, INPUT = INPUT, clean = clean))
				os.system('sleep 20')


MACS2('wt')
MACS2('ko')

	


