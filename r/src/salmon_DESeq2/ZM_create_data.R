
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")
data_directory = file.path(getwd(), glue('data/organised/{treatment}/output_salmon'))

times = c(0, 16, 24, 36, 48)


dds <- create_dds('ZM', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "ZM_data.RData")))
        