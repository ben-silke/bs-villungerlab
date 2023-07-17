
library("vsn")
library("genefilter")
setwd("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

times = c(0, 8, 12, 16, 24, 48)
treatment <- "Etop"
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))


dds_Etop_r1to6 <- create_dds('Etop', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds_Etop_r1to6, file = glue('r/data/', "Etop_r1to6_data.RData"))

res <- results(dds_Etop_r1to6)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results/Etop_r1to6_data.csv")
