
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")
treatment <- "Etop"
data_directory = file.path('/Volumes/bs_external/villunger', glue('organised/{treatment}/output_salmon'))

times = c(0, 8, 12, 16, 24, 48)


dds <- create_dds('Etop', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds, file = glue('r/data/', "Etop_r1to6_data.RData"))

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results/Etop_r1to6_data.csv")

        