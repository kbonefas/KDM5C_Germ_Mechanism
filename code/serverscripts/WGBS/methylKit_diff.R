# 24.08.24 DNA methylation analysis of bismark bam files

#read in the bam files
bamfiles <- list.files(path = "../data/processed/BAM_dedup/", pattern = "*.sorted.dedup.bam",
           full.names = TRUE)

print(bamfiles)

# my.methRaw=processBismarkAln(location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
#                              package = "methylKit"),
# 		                         sample.id="test1", assembly="mm10", 
# 		                         read.context="CpG", save.folder=getwd())