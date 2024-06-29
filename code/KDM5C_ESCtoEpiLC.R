
WesternQuant <- function(csvfile, expprotein, housekeep){
    data <- read.csv(file = csvfile, sep = ",", header = T)
    colnames(data) <-c("Channel","Signal","Genotype","Time","Area","Mean","Min","Max")
    head(data)

    #subtract background from signal for normalized value
    data_green_sig <- subset(data, Channel == "green" & Signal == expprotein)
    data_green_bg <- subset(data, Channel == "green" & Signal == "background")
    data_green <- data.frame(data_green_sig, Normalized = data_green_sig$Mean - data_green_bg$Mean)

    data_red_sig <- subset(data, Channel == "red" & Signal == housekeep)
    data_red_bg <- subset(data, Channel == "red" & Signal == "background")
    data_red <- data.frame(data_red_sig, Normalized = data_red_sig$Mean - data_red_bg$Mean)

    Loading <- data_green$Normalized/data_red$Normalized
	data_green[data_green == "5cKO"] <- "5CKO"
    df <- data.frame(data_green, Loading, GenoTime = paste(data_green$Genotype, data_green$Time) )

    #set the order
    df$GenoTime <- factor(df$GenoTime, levels = c("WT 0", "WT 24", "WT 48","5CKO 0", "5CKO 24", "5CKO 48"))

    library(ggpubr)
    #check the dataframe looks correct
    print(df)
	nrow(df)

    my_comparisons <- list(c("WT 0", "WT 24"), c("WT 24", "WT 48"))
    pl <- ggboxplot(df, x = 'GenoTime', y = 'Loading', color = "black", fill="GenoTime", 
        title = "KDM5C protein expression", 
        #subtitle = "ESC to EpiLC Differentiation",
        add = "dotplot", add.params = list(size = 1.5),
        xlab = " ", ylab = "Normalized KDM5C/DAXX Intensity [au]",
		ylim = c(0, 2),
        palette = c("WT 0" = "#f2b722ff", "WT 24" = "#ef8e3eff", "WT 48" = "#b26c93ff","5CKO 0" = "#f2b722ff", "5CKO 24" = "#ef8e3eff", "5CKO 48" = "#b26c93ff")) +
        stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif", size = 6) +
        rremove("legend") +
		scale_x_discrete(labels = c("0\nWT", "24\nWT", "48\nWT", "0\n5CKO", "24\n5CKO", "48\n5CKO"))
    

    ggsave(snakemake@output[[1]], plot = pl, width = 4, height = 4)

}

WesternQuant(snakemake@input[[1]], snakemake@params[["expprotein"]], snakemake@params[["housekeep"]])
