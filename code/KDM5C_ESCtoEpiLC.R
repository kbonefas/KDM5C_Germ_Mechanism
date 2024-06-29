#24.06.24 KDM5C western ESC to EpiLC expression

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


#24.06.29 KDM5C RNA ESC to EpiLC expression

Kdm5c_qpcr <- read.csv(snakemake@input[[2]], sep =",", header = T)
colnames(Kdm5c_qpcr) <- c("Sample", "Genotype", "Timepoint", "Primer","CT", "Replicate")
head(Kdm5c_qpcr)


#get the average CT value for each technical replicate
avg_CT <- function(df, samp, primey){
  reps <- subset(df, Sample == samp & Primer == primey)
  avg = (reps[1,5]+reps[2,5])/2
  return(avg)
}

#difference between CT values for replicates to see if there are any major deviations
diff_CT <- function(df, samp, primey){
  reps <- subset(df, Sample == samp & Primer == primey)
  diff = abs(reps[1,5]-reps[2,5])
  return(diff)
}

#get the information that goes with the average (sample, timepoint, primer etc)
info_CT <- function(df, samp, primey){
  reps <- subset(df, Sample == samp & Primer == primey)
  rowsamp = reps[1,1:4]
  return(rowsamp)
}

#checking to see if it works
info_CT(Kdm5c_qpcr, "0-2a", "TBP")


#list of all primers and samples
primerlist <- unique(Kdm5c_qpcr$Primer)
samplelist <- unique(Kdm5c_qpcr$Sample)

cnt <- 1

#go through each primer and sample and average the replicates together. Store them in a list and dataframe
CT_avg <- c()
outsamp <- data.frame()
CT_diff <- c()

for (p in primerlist){
  for (s in samplelist){
    CT_avg[cnt] <- avg_CT(Kdm5c_qpcr, s, p)
    CT_diff[cnt] <- diff_CT(Kdm5c_qpcr, s, p)
    cnt = cnt + 1
    
    rowsamp <- info_CT(Kdm5c_qpcr, s, p)
    outsamp = rbind(outsamp, rowsamp)
  }
}

Kdm5c_avgs <- data.frame(outsamp, CT_avg, CT_diff)
head(Kdm5c_avgs)

#remove the controls 
Kdm5c <- subset(Kdm5c_avgs, !is.na(Timepoint))

#####2^-dct for Kdm5c against TBP
#SEM
std <- function(x) sd(x)/sqrt(4)
#2^-dCT
twotodct <- function(a,b) 2^(-(a-b))
#SEM <- (std(res.1$Ct,res.2$Ct))
twodct_over_TBP <- function(df, samp, primey){
  sampy <- subset(df, Sample == samp & Primer == primey)
  TBP <- subset(df, Sample == samp & Primer == "TBP")
  twonegdct = twotodct(sampy[1,5],TBP[1,5])
  return(twonegdct)
}


samplelist_2 <- unique(Kdm5c$Sample)

#store results in empty list and dataframe
twonegCT_TBP <- c()
outdCT <- data.frame()

count = 1

for (s in samplelist_2){
  twonegCT_TBP[count] <- twodct_over_TBP(Kdm5c, s, "KDM5C")
  count = count + 1
  
  rowsamp <- info_CT(Kdm5c, s, "KDM5C")
  outdCT = rbind(outdCT, rowsamp)
}


tail(outdCT)


Kdm5c_dCT_TBP <- data.frame(outdCT, twonegCT_TBP) 

Kdm5c_dCT_TBP[Kdm5c_dCT_TBP == "5cKO"] <- "5CKO"
Kdm5c_dCT_TBP$GenoTime <- paste(Kdm5c_dCT_TBP$Genotype, Kdm5c_dCT_TBP$Timepoint)
tail(Kdm5c_dCT_TBP)
#############Make boxplot
#set the order
Kdm5c_dCT_TBP$GenoTime <- factor(Kdm5c_dCT_TBP$GenoTime, levels = c("WT 0", "WT 24", "WT 48","5CKO 0", "5CKO 24", "5CKO 48"))


my_comparisons <- list(c("WT 0", "WT 24"), c("WT 24", "WT 48"))
pl <- ggboxplot(Kdm5c_dCT_TBP, x = 'GenoTime', y = 'twonegCT_TBP', color = "black", fill="GenoTime", 
    title = "Kdm5c mRNA expression", 
    #subtitle = "ESC to EpiLC Differentiation",
    add = "dotplot", add.params = list(size = 1.5),
    xlab = " ", ylab = "2^-dCT (normalized to TBP)",
	ylim = c(0, 0.2),
    palette = c("WT 0" = "#f2b722ff", "WT 24" = "#ef8e3eff", "WT 48" = "#b26c93ff","5CKO 0" = "#f2b722ff", "5CKO 24" = "#ef8e3eff", "5CKO 48" = "#b26c93ff")) +
    stat_compare_means(comparisons = my_comparisons, method="t.test", label = "p.signif", size = 6) +
    rremove("legend") +
	scale_x_discrete(labels = c("0\nWT", "24\nWT", "48\nWT", "0\n5CKO", "24\n5CKO", "48\n5CKO"))
    

ggsave(snakemake@output[[2]], plot = pl, width = 4, height = 4)



#WT samples only
# Kdm5c_dCT_TBP.wt <- subset(Kdm5c_dCT_TBP, Genotype == "WT")
# 
# pl.wt <- ggplot(Kdm5c_dCT_TBP.wt, aes(x=GenoTime, y=twonegCT_TBP, fill=GenoTime)) +
#   geom_boxplot(position=position_dodge(0.8))+
#   geom_dotplot(binaxis='y', stackdir='center', 
#                position=position_dodge(0.8), dotsize=1) +
#   # ylim(c(-0.1, 1.5))+
#   ggtitle("Kdm5C ESC to EpiLC mRNA Expression") +
#   xlab(" ") + ylab("2^-dCT (normalized to TBP)") + 
#   scale_fill_manual("Genotype", values = c("WT 0" = "#f2b722ff", "WT 24" = "#ef8e3eff", "WT 48" = "#b26c93ff","5cKO 0" = "#f2b722ff", "5cKO 24" = "#ef8e3eff", "5cKO 48" = "#b26c93ff"))+
#   theme(text = element_text(size = 20))   
# 
# 
# pdf(file = "KDM5C_mRNA_ESCtoEpiLC_WTonly.pdf",   # The directory you want to save the file in
#     width = 9, # The width of the plot in inches
#     height = 6) # The height of the plot in inches
# pl.wt
# dev.off()


##statistics
#Welch's ttest
Kdm5c_dCT_TBP.wt.0 <- subset(Kdm5c_dCT_TBP, Timepoint == "0" & Genotype == "WT")
Kdm5c_dCT_TBP.wt.24 <- subset(Kdm5c_dCT_TBP, Timepoint == "24" & Genotype == "WT")
Kdm5c_dCT_TBP.wt.48 <- subset(Kdm5c_dCT_TBP, Timepoint == "48" & Genotype == "WT")

Kdm5c_024 <- t.test(Kdm5c_dCT_TBP.wt.0$twonegCT_TBP,Kdm5c_dCT_TBP.wt.24$twonegCT_TBP)
Kdm5c_2448 <- t.test(Kdm5c_dCT_TBP.wt.24$twonegCT_TBP,Kdm5c_dCT_TBP.wt.48$twonegCT_TBP)
Kdm5c_048 <- t.test(Kdm5c_dCT_TBP.wt.0$twonegCT_TBP,Kdm5c_dCT_TBP.wt.48$twonegCT_TBP)


