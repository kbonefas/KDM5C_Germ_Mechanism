#2026.07.20 - make the sample sheet for consensus peaks

#1) #make a dataframe with the sample infomration

names <- c("KDM5C_rep1", "KDM5C_rep2")

sample_sheet <- data.frame(SampleID = names, 
    Tissue = rep("ESC", 2),
    Factor = rep("KDM5C", 2),
    Condition = rep("KDM5C", 2),
    Replicate = c("1", "2"),
    bamReads = c("../../data/kdm5c/bam/sort_Agarwal_KDM5C_rep1.bam", "../../data/kdm5c/bam/sort_Agarwal_KDM5C_rep2.bam"),
    bamControl = rep("../../data/kdm5c/bam/sort_Agarwal_Input_merged.bam", 2),
    Peaks = c("../../data/kdm5c/macs2/Agarwal_KDM5C_rep1_filtered_peaks.narrowPeak", "../../data/kdm5c/macs2/Agarwal_KDM5C_rep2_filtered_peaks.narrowPeak"),
    PeakCaller = rep("macs", 2)
)

print(head(sample_sheet))

write.table(sample_sheet, "../../data/kdm5c/macs2/Agarwal_ESC_sample_sheet.csv", quote = FALSE, sep = ",", row.names = FALSE)

## 26.07.20 - get consensus peaks from MACS3 ChIPseq bed files - Snakemake
#made script more robust
# to run from the command line using the r environment:
    # conda activate r_env
    # R -e "source('ChIPseq_consensus_peaks.R')"


library(AnnotationDbi)
library(GO.db)
library(DiffBind)


#function to get the consensus peaks
#name is how you differentiate the datasets, for example EpiLC or PNC"
#samplesheet - the sample sheet name
print("working directory:")
print(getwd())

peaks <- dba(sampleSheet=sample_sheet)
print(head(peaks))

# Step 1: Generate the consensus peaksets from your conditions
peaks <- dba.peakset(peaks, consensus = DBA_CONDITION, minOverlap = 1)

# Step 2: Use dba() to mask out everything except the consensus sets
consensus_only_object <- dba(peaks, mask = peaks$masks$Consensus)

# Step 3: Retrieve the clean consensus peak intervals as a GRanges object
consensus_peaks <- dba.peakset(consensus_only_object, bRetrieve = TRUE)
consensus_peaks <- as.data.frame(consensus_peaks)
print("consensus_peaks:")
print(head(consensus_peaks))

#just the bedfile coordinates
consensus_peaks_df <- consensus_peaks[c(1:3)]
print(head(consensus_peaks_df))
write.table(consensus_peaks_df, file="../../data/kdm5c/macs2/Agarwal_ESC_KDM5C_consensus_peaks.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


