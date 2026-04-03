### 2026.03.25 - CpG Germ

# BiocManager::install(c("GenomicRanges", "GenomicFeatures", 
#                        "BSgenome.Mmusculus.UCSC.mm10", 
#                        "Biostrings", "EnrichedHeatmap"))


# BiocManager::install("Biostrings")


library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(EnrichedHeatmap)

# Load genome
genome <- BSgenome.Mmusculus.UCSC.mm10

getwd()

# Load gene annotation (GTF)
txdb <- makeTxDbFromGFF("../KDM5C_Resubmission/gencode.vM23.annotation.gtf")

# Get promoters (TSS ± 2kb)
promoters <- promoters(genes(txdb), upstream=500, downstream=500)

# Tile promoters into small bins (e.g., 50 bp)
BINSIZE = 3
bins <- unlist(tile(promoters, width=BINSIZE))

# Extract sequences
seqs <- getSeq(genome, bins)

# Function to compute CpG observed/expected
cpg_oe <- function(seq) {
  seq <- toupper(as.character(seq))
  c_count <- stringr::str_count(seq, "C")
  g_count <- stringr::str_count(seq, "G")
  cg_count <- stringr::str_count(seq, "CG")
  len <- nchar(seq)
  
  if (c_count == 0 || g_count == 0) return(0)
  (cg_count * len) / (c_count * g_count)
}

# Compute CpG values
library(stringr)
cpg_values <- sapply(seqs, cpg_oe)
print(head(cpg_values))

# Convert to matrix (rows = promoters, cols = bins)
n_bins <- 4000 / BINSIZE  # 4kb window / 50bp bins
mat <- matrix(cpg_values, ncol=n_bins, byrow=TRUE)
print(head(mat))

# Plot profile
pdf(snakemake@output[[1]], width = 4, height = 4)

# Make plots
plot(colMeans(mat), type="l", lwd=2,
     xlab="Position relative to TSS",
     ylab="CpG observed/expected",
     main="CpG content across promoters")

abline(v=n_bins/2, col="red", lty=2)  # TSS position

dev.off()

