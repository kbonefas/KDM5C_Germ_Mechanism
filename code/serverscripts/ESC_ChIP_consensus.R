#2026.07.20 - make the sample sheet for consensus peaks

#1) #make a dataframe with the sample infomration

names <- c("KDM5C_rep1", "KDM5C_rep2")



sample_sheet <- data.frame(SampleID = names, 
    Tissue = rep("ESC", 2),
    Factor = rep("KDM5C", 2),
    Condition = rep("KDM5C", 2),
    Replicate = c("1", "2"),
    bamReads = c(),
    bamControl = rep(Controls, 2),
    Peaks = snakemake@input[["peaks"]],
    PeakCaller = rep("macs", 2)
)

print(head(sample_sheet))

write.table(sample_sheet, "Agarwahl_ESC_sample_sheet.csv", quote = FALSE, sep = ",", row.names = FALSE)

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

consensus_peaks <- dba.peakset(peaks, consensus = DBA_CONDITION, minOverlap = 1, bRetrieve = TRUE)
consensus_peaks <- as.data.frame(consensus_peaks)
print("consensus_peaks:")
print(head(consensus_peaks))

#just the bedfile coordinates
consensus_peaks_df <- consensus_peaks[c(1:3)]
print(head(consensus_peaks_df))
write.table(consensus_peaks_df, file="Agarwahl_ESC_KDM5C_consensus_peaks.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


