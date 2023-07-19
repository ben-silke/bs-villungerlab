

library("vsn")
library("genefilter")
library(dplyr)

setwd("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

short_times <- c(8,12,16,24,48)
long_times <- c(16,20,24,36,48)

# Create a named list
times <- list(
  'ZM' = long_times,
  'Noc' = long_times,
  'DHCB' = long_times,
  'Nutl' = short_times,
  'Etop' = short_times
)

treatments <- c("ZM", "DHCB", "Etop", "Noc", "Nutl")
path <- '../../../Volumes/bs_external/lab_villunger'

data_frame <- data.frame(treatments = treatments)

dds_es = list()

for (treatment in treatments) {
    data_directory = file.path(path, glue('organised/{treatment}/output_salmon'))
    dds <- create_dds(treatment, data_directory, times[[treatment]], "salmon_quant", 1:6)
    dds_es <- append(dds)
}
data_frame$ddses <- dds_es


