
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")
data_directory = file.path('/Volumes/bs_external/villunger', glue('organised/{treatment}/output_salmon'))

times = c(0, 16, 20, 24, 36, 48)
treatment <- "ZM"


dds <- create_dds('ZM', data_directory, times, "salmon_quant", 1:3)
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "ZM_r1to3_data.RData")))

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results/ZM_r1to3_data.csv")

        