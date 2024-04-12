#############################################################################################
#                                                                                           #
# 2022 10 20 For generating a basic MA plot using ggpubr                                    # 
#                                                                                           #
# Inputs: DESeq2 results file (.csv), padj cut off                                          #
# Libraries: ggpubr                                                                         #
# Outputs: pdf of                                                                      #
#                                                                                           #
#############################################################################################

genmaplot <- function(DESEQ){
  ggmaplot(DESEQ, 
           fdr = 0.1, fc = 0, size = 0.6,
           palette = c("#f2b722ff", "#1465AC", "#696969"),
           genenames = as.vector(DESEQ$symbol),
           legend = "top", top = 10,
           font.label = c("bold", 11), label.rectangle = TRUE,
           font.legend = "bold",
           # ylim = c(-3.5, 3.5),
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())
}