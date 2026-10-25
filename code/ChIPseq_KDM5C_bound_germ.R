#2026.10.25 KDM5C binding at germline genes across all datasets (pESC, EpiLC, NPC)


samples <- c("ESC","EpiLC", "PNC") #make sure order matches snakefile

germ <- read.csv(snakemake@input[[3]], sep = ",")
ESC_bound <- read.csv(snakemake@input[[1]], sep = ",")
EpiLC_bound <- read.csv(snakemake@input[[2]], sep = ",")


germ$Kdm5c_binding <- ifelse(germ$ENSEMBL %in% ESC_bound$ENSEMBL & germ$ENSEMBL %in% EpiLC_bound$ENSEMBL, "pESC&EpiLC", 
    ifelse(germ$ENSEMBL %in% ESC_bound$ENSEMBL & !(germ$ENSEMBL %in% EpiLC_bound$ENSEMBL), "pESC_only", 
        ifelse(germ$ENSEMBL %in% EpiLC_bound$ENSEMBL & !(germ$ENSEMBL %in% ESC_bound$ENSEMBL), "EpiLC_only", "Unbound")))

print(head(germ))

write.table(germ, snakemake@output[[1]], sep = ",", row.names = FALSE)