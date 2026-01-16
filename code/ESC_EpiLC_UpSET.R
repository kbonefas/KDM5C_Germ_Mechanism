#2023.01.13 - Upset plot germline genes 
germ <- read.csv(snakemake@input[[1]])
print("germline genes")
print(head(germ))

DEGs <- list()
#read in the DEGs
samples <- c("ESC", "48NO", "48VA", "96NO", "96VA")
for(i in 1:length(samples)){
    df <- read.csv(snakemake@input[[1+i]], sep = ",", row.names = 1)
    print(head(df))

    DEGs[[i]] <- row.names(df)
}

names(DEGs) <- samples


#### plot simple upset
library("UpSetR")

modifiedupset <- function(samplelist){
	upset(fromList(samplelist), order.by = "freq",  sets.x.label = "# Germline DEGs", mainbar.y.label = "# of Overlapping Germline DEGs", empty.intersections = "on")
}

pdf(file = snakemake@output[[1]], width = 8, height = 5)

upset(fromList(DEGs), order.by = "freq",  sets.x.label = "# germline DEGs", mainbar.y.label = "# in group", text.scale = 2, sets = samples, mb.ratio = c(0.55, 0.45), keep.order = TRUE)

dev.off()

#vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
#need a cagetory row in my upset list??? that has the category KDM5C bound
    #are you even able to do it in list format
# apparently you can do scatter plots with UpsetR
# Comparing two sets in a whole: https://github.com/hms-dbmi/UpSetR?tab=readme-ov-file



#plotting
library("ggplot2")







# library(ComplexUpset)
# movies = as.data.frame(ggplot2movies::movies)
# genres = colnames(movies)[18:24]

# # for simplicity of examples, only use the complete data points
# movies[movies$mpaa == '', 'mpaa'] = NA
# movies = na.omit(movies)


# upset(
#     movies,
#     genres,
#     base_annotations=list(
#         'Intersection size'=intersection_size(
#             counts=FALSE,
#             mapping=aes(fill=mpaa)
#         )
#     ),
#     width_ratio=0.1
# )

# library(UpSetR)
# test <- as.data.frame(matrix(rnorm(1:12), nrow = 4))
# colnames(test) <- c("A", "B", "C")
# test_up <- test
# test_up[test >= 0] <- 1
# test_up[test_up != 1] <- 0
# test_up$Gene <- paste0("Gene", 1:nrow(test), "_UP")
# test_down <- test
# test_down[test < 0] <- 1
# test_down[test_down != 1] <- 0
# test_down$Gene <- paste0("Gene", 1:nrow(test), "_DOWN")
# test <- rbind(test_up, test_down)


# upordown <- function(row, direction) {
#   gene <- row["Gene"]
#   if (grepl(x = gene, pattern = direction)) {
#     newData <- T
#   }
#   else {
#     newData <- F
#   }
# }


# metadata <- data.frame(
#   c("A", "B", "C"),
#   as.numeric(apply(test_up[, 1:3], 2, sum))
# )
# colnames(metadata) <- c(
#   "sets",
#   "NumberUP"
# )


# upset(test,
#       sets = c("A", "B", "C"), set.metadata = list(
#         data = metadata,
#         plots = list(
#           list(type = "hist", 
#                column = "NumberUP", 
#                assign = 20, # defines width of the meta-data histogram
#                colors = "red")
#         )
#       ),
#       queries = list(list(query = upordown, params = list("_UP"), color = "red", active = TRUE))
# )