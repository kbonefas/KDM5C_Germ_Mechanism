#2026.07.21 Filtering the narrowpeak files for high confidence LHX2 peaks

# #example rule
# rule LHX2_confidence:
#     input:
#         "data/processed/MACS3/{LHX2sample}_narrow_peaks.narrowPeak"
#     params:
#         cutoff = 5
#     output:
#         "data/processed/MACS3/{LHX2sample}_narrow_peaks_filtered.narrowPeak"
#     script:
#         "code/utilities/MACS3_narrowPeak_filter.R"


#read in narrowpeak file
peaks <- c("../../data/kdm5c/macs2/Agarwal_KDM5C_rep1_peaks.narrowPeak", 
            "../../data/kdm5c/macs2/Agarwal_KDM5C_rep2_peaks.narrowPeak")
for (i in 1:length(peaks)){
    narrow <- read.csv(peaks[i], sep = "\t", header = FALSE)
colnames(narrow) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
#signalValue = Measurement of overall (usually, average) enrichment for the region. In MACS3, we use the fold-enrichment values here. Please note that, same as the columns 8 and 9, this value corresponds to the peak summit, where the enrichment is the highest in the peak region.

# print(head(narrow))



#subset for loci that have a signalValue greater than the cut off
narrow_sub <- subset(narrow, narrow$qValue >= 10)

#sort the peak file by lowest to highest qValue
narrow_order <- narrow_sub[order(narrow_sub$qValue),]

print("lowest")
print(head(narrow_order, 10))

#sort the peak file by highest to lowest qValue
narrow_order <- narrow_sub[order(-narrow_sub$qValue),]
print("highest")
print(head(narrow_order, 10))


write.table(narrow_sub, paste0("../../data/kdm5c/macs2/Agarwal_KDM5C_rep", i, "_filtered_peaks.narrowPeak"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
