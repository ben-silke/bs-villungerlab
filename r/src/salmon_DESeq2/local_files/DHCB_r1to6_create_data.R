
library("vsn")
library("genefilter")
setwd("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

times = c(0, 16, 20, 24, 36, 48)
treatment <- "DHCB"
data_directory = file.path('/nobackup/lab_villunger/bsilke', glue('organised/{treatment}/output_salmon'))

dds_DHCB_r1to6 <- create_dds('DHCB', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds_DHCB_r1to6, file = glue('r/data/', "DHCB_r1to6_data.RData"))

res <- results(dds_DHCB_r1to6)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results/DHCB_r1to6_data.csv")
