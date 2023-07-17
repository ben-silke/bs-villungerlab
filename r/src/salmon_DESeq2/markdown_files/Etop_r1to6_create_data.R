
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

res_Etop_r1to6 <- results(dds_Etop_r1to6)
resOrdered_Etop_r1to6 <- res_Etop_r1to6[order(res_Etop_r1to6$padj),]
resOrdered_Etop_r1to6 <- add_annotations_to_results(resOrdered_Etop_r1to6)

head(resOrdered_Etop_r1to6)

resOrderedDF_Etop_r1to6 <- as.data.frame(resOrdered_Etop_r1to6)
write.csv(resOrderedDF_Etop_r1to6, file = "results/Etop_r1to6_data.csv")
