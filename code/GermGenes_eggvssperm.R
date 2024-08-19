### 24.06.11 Define male (sperm) and female (egg) germline genes and compare XX vs XY EpiLCs

# Read in germline gene list 
germ <- read.csv(snakemake@input[[1]], sep = ",")

# Sankey diagram of XX vs XY germline genes

library(ggsankey)
library(ggplot2)
library(dplyr)


sankeydf <- data.frame(ENSEMBL = germ$ENSEMBL, Total = rep("germline genes", nrow(germ)), XX = ifelse(germ$sexBias == "XX", "egg-biased", "no"), XY = ifelse(germ$sexBias == "XX", NA, ifelse(germ$sexBias == "XY", "sperm-biased", "unbiased")))


# sankey_long <- data.frame(row.names = germ$ENSEMBL, )

#XY = ifelse(germ$sexBias == "XY", "sperm-biased", "unbiased"))
#sankeydf <- data.frame(ENSEMBL = germ$ENSEMBL, Total = rep("germ genes", nrow(germ)), XX = ifelse(germ$sexBias == "XX", "egg-biased", "no"), XY = ifelse(germ$sexBias == "XX", "egg-biased", ifelse(germ$sexBias == "XY", "sperm-biased", ifelse(germ$sexBias == "unbiased", "unbiased", )))


row.names(sankeydf) <- sankeydf$ENSEMBL
sankeydf <- subset(sankeydf, select = -ENSEMBL)

print(head(sankeydf))

##2) make the dataframe for plotting, including how many observations there are for each category
df <- sankeydf %>%
  make_long(Total, XX, XY)

print(head(df))

# count how many there are in each group and merge back with the long df
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
print("plotting df")
print(head(df2))
print(unique(df2$node))
print(unique(df2$next_node))

df2$node <- factor(df2$node, levels = c("germline genes", "no", "unbiased", "sperm-biased", "egg-biased" ))
df2$next_node <- factor(df2$next_node, levels = c("no", "unbiased", "sperm-biased", "egg-biased" ))

# plot the dataframe

pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = paste0(node," n=", n))) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.2) +
  scale_fill_manual(values = c("#5fa35f", "#9e3c8e", "mediumpurple1", "dodgerblue1", "hotpink1"), drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Egg vs Sperm-biased Mouse Germline-enriched Genes")


ggsave(snakemake@output[[1]], plot = pl, width = 12, height = 3)



### Make bar graphs of the sperm vs egg biased DEGs
##1) Read in the EpiLC XX and XY germline DEGs
source("code/utilities/parameters.R")
samples <- c("XY5cKO", "XX5cHET", "XX5cKO")
DEGs <- list()
for (s in 1:length(samples)){
	#read in the dataframe
	df <- read.csv(snakemake@input[[s+1]], sep = ",")
	
	gdf <- subset(df, ENSEMBL %in% germ$ENSEMBL & log2FoldChange > l2fcco_ESCEpiLC)
	germGENES <- gdf$ENSEMBL
	#subset for the number of germline genes
	DEGs[[s]] <- germGENES
}
names(DEGs) <- samples


##1) Separate out number of genes in each bias
#DEGs is list of germline
egg_germ <- subset(germ, sexBias == "XX")
sperm_germ <- subset(germ, sexBias == "XY")
unbiased_germ <- subset(germ, sexBias == "unbiased")

#create dataframe for plotting
qlotdf <- data.frame()

for (s in samples){
	#read in the dataframe
	a <- DEGs[[s]]
	egg_bias <- length(a[a %in% egg_germ$ENSEMBL])
	sperm_bias <- length(a[a %in% sperm_germ$ENSEMBL])
	Unbiased <- length(a[a %in% unbiased_germ$ENSEMBL])
	df <- data.frame(Genotype = rep(s, 3), SexBias = c("Egg-biased", "Sperm-biased", "Unbiased"), Count = c(egg_bias, sperm_bias, Unbiased))
	qlotdf <- rbind(qlotdf, df)
}

print(qlotdf)
#set the plotting order
qlotdf$SexBias <- factor(qlotdf$SexBias, levels = c("Egg-biased", "Unbiased", "Sperm-biased"))

samples <- c("XY5cKO", "XX5cHET", "XX5cKO")

qlotdf$Genotype[qlotdf$Genotype == 'XY5cKO'] <- 'XY 5CKO'
qlotdf$Genotype[qlotdf$Genotype == 'XX5cHET'] <- 'XX 5CHET'
qlotdf$Genotype[qlotdf$Genotype == 'XX5cKO'] <- 'XX 5CKO'


qlotdf$Genotype <- factor(qlotdf$Genotype, levels = c("XY 5CKO", "XX 5CHET", "XX 5CKO"))


##2) plot the data
library("ggpubr")
q <- ggbarplot(qlotdf, "Genotype", "Count",
  fill = "SexBias", color = "SexBias", palette = c("Egg-biased" = "hotpink1", "Sperm-biased" = "dodgerblue1", "Unbiased" = "mediumpurple1" ),
  label = TRUE, lab.col = "black", lab.vjust = 2, lab.hjust = 0.5, xlab = " ", ylab = "# of DEGs", orientation = "vert") 

q <- ggpar(q, legend = "bottom", legend.title = "Sex-Biased Expression in Germline", main = "All Germline DEGs")


### 23.11.13 Make a bar graph with the proportion of UNIQUE germline DEGs in each sample that are XX vs XY vs unbiased
XXallDEGs <- c(DEGs[["XX5cHET"]], DEGs[["XX5cKO"]])
XX_UNIQUE <- setdiff(XXallDEGs, DEGs[["XY5cKO"]])
XY_UNIQUE <- setdiff(DEGs[["XY5cKO"]], XXallDEGs)
uniquelist <- list(XX_UNIQUE, XY_UNIQUE)

sample2 <- c("XX_only", "XY_only")
names(uniquelist) <- sample2
qlotdf2 <- data.frame()

for (s in sample2){
	#read in the dataframe
	a <- uniquelist[[s]]
	Egg_bias <- length(a[a %in% egg_germ$ENSEMBL])
	Sperm_bias <- length(a[a %in% sperm_germ$ENSEMBL])
	Unbiased <- length(a[a %in% unbiased_germ$ENSEMBL])
	df <- data.frame(Sex = rep(s, 3), SexBias = c("Egg-biased", "Sperm-biased", "Unbiased"), Count = c(Egg_bias, Sperm_bias, Unbiased))
	qlotdf2 <- rbind(qlotdf2, df)
}

print(qlotdf2)
#set the plotting order
qlotdf2$Sex <- factor(qlotdf2$Sex, levels = sample2)
qlotdf2$SexBias <- factor(qlotdf2$SexBias, levels = c("Egg-biased", "Sperm-biased", "Unbiased"))

##2) plot the data
library("ggpubr")
t <- ggbarplot(qlotdf2, "Sex", "Count",
  fill = "SexBias", color = "SexBias", palette = c("Egg-biased" = "hotpink1", "Sperm-biased" = "dodgerblue1", "Unbiased" = "mediumpurple1" ),
  label = TRUE, lab.col = "black", lab.vjust = 2, lab.hjust = 0.5, xlab = " ", ylab = "# of DEGs", orientation = "vert") + scale_x_discrete(labels = c("XX Unique", "XY Unique"))

t <- ggpar(t, legend = "bottom", main = "Sex-specific DEGs")

library("gridExtra")
eggsperm <- list(q, t)
ggsave(snakemake@output[[2]], plot = grid.arrange(grobs = eggsperm, nrow = 1, ncol = 2), width = 7, height = 3.5)
