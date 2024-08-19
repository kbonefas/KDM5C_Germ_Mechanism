#### 23.10.30 fastq alignment using BOWTIE1 to mm10
# first activate the blueberry environment in conda so python is on version 3
	#you might have to deactivate from the base enrvironment first
import os
import glob
import time

WigToBigwig = "/home/saurabha/UTILITIES/UCSC/wigToBigWig"


#algining to mm10
#folders with alginment information
GENOME_NAME = "mm10"
GENOME_FOLDER = "/nfs/value/siwase2/GENOMES/MM10"
GENOME_FASTA = GENOME_FOLDER + "/{GENOME_NAME}.fa".format(GENOME_NAME = GENOME_NAME)
GENOME_FaFai = GENOME_FOLDER + "/{GENOME_NAME}.fa.fai".format(GENOME_NAME = GENOME_NAME)
GENOME_INDEX = "/nfs/value/siwase2/GENOMES/MM10/mm10_bowtie1"

## MAPPING ##
#where bowtie is located
BOWTIE="/home/saurabha/MAPPING_SOFTWARES/BOWTIE_1/bowtie-1.1.2/bowtie"



#samples - fastq.gz files in the data/raw folder
samples =  glob.glob("../data/raw/*.fastq.gz")
for s in samples:
	#name of sample without extensions, adding in genome
	name = s.replace("../data/raw/", "")
	name = name.replace("fastq.gz", GENOME_NAME)
	
	#run the bowtie commands to align
	os.system("gzip -cd {sample} | {BOWTIE} -v 2 -m 1 --best --strata -S --time -p 16 {GENOME_INDEX} - | samtools view -F 0x4 -ub -t {GENOME_FaFai} |samtools sort -m 192G -@ 10 - -o ../data/processed/BAM/Sorted_{name}.bam".format(sample = s, BOWTIE = BOWTIE, GENOME_INDEX = GENOME_INDEX, GENOME_FaFai = GENOME_FaFai, name = name))
	
	#remove the duplicates
	os.system("samtools view -F 0x4 -ub ../data/processed/BAM/Sorted_{name}.bam | samtools rmdup -s --output-fmt-option nthreads=8 - ../data/processed/BAM/{name}.nodup.bam".format(name = name))
	
	#index the bam file
	os.system("samtools index ../data/processed/BAM/{name}.nodup.bam".format(name = name))
	time.sleep(10)

	## READ COUNTING ##
	#count the number of reads in each alignment
	os.system("echo {name} >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(name = name, GENOME_NAME = GENOME_NAME))
	os.system("echo Mapped >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(GENOME_NAME = GENOME_NAME))
	awkcommand = "awk -F '\t' '{sum+=$3;} END{print sum;}'"
	os.system("samtools idxstats ../data/processed/BAM/{name}.nodup.bam | {awkcommand} >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(name = name, awkcommand=awkcommand, GENOME_NAME = GENOME_NAME))
	
	#count the reads aligned to ribosomes and mitochondria
	os.system ("echo Ribosomal >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(GENOME_NAME = GENOME_NAME))
	os.system("samtools view -c ../data/processed/BAM/{name}.nodup.bam chr17:39978942-39986774 >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(name = name, GENOME_NAME = GENOME_NAME))
	os.system ("echo Mitochondrial >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(GENOME_NAME = GENOME_NAME))
	os.system("samtools view -c ../data/processed/BAM/{name}.nodup.bam chrM >> ../data/processed/BAM/READ_COUNTS_{GENOME_NAME}.txt".format(name = name, GENOME_NAME = GENOME_NAME))

