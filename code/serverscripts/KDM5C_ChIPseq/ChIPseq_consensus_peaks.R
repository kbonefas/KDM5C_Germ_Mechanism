## 24.04.26 - get consensus peaks for filtered ChIPseq bed files 
### 24.05.14 update - only do for PNCs, as EpiLC data is only one replicate
#made script more robust
# to run from the command line using the r environment:
    # conda activate r_env
    # R -e "source('ChIPseq_consensus_peaks.R')"


library(AnnotationDbi)
library(GO.db)
library(DiffBind)

print("current working directory:")
print(getwd())

#where the output files are going to go
outpath <- '../data/processed/'

#function to get the consensus peaks
#name is how you differentiate the datasets, for example EpiLC or PNC"
#samplesheet - the sample sheet name

consensus <- function(name, samplesheet){
    #read in the sample sheet file, using the file names of the new filtered narrowpeaks files
    KDM5C_peaks <- read.csv(samplesheet)
    #get the consensus peaks for WT and KO separately

    #get the conditions (genotypes) in the sample sheet
    #for both WT and KO, peaks can be in only one replicate, since we're removing KO peaks later
    genos <- unique(KDM5C_peaks$Condition)

    for (i in genos){
        geno_samples <- subset(KDM5C_peaks, Condition == i)
        peaks_KDM5C_geno <- dba(sampleSheet= geno_samples)
        print(head(peaks_KDM5C_geno))


        consensus_peaks <- dba.peakset(peaks_KDM5C_geno, consensus = DBA_CONDITION, minOverlap = 1, bRetrieve = TRUE)
        consensus_peaks <- as.data.frame(consensus_peaks)
        consensus_peaks_df <- consensus_peaks[c(1:3)]
        #consensus_peaks_df_WT$seqnames <- as.character(as.factor(consensus_peaks_df_WT$seqnames))
        print(head(consensus_peaks_df))
        write.table(consensus_peaks_df, file= paste0(outpath,'consensuspeaks_mm10_', name,'_',i,'.bed'), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

    }

}

consensus('PNC', 'samplesheet_mm10_KDM5C_PNC.csv')
