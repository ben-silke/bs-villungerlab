
getwd()
source("r/src/utils.R")
source("r/src/salmon_pca_analysis/pca_analysis_util.R")

treatment <- "ZM"
salmon_data_directory = file.path(getwd(), glue('data/organised/{treatment}/output_salmon'))
times = c(0, 16, 20, 24, 36, 48)

print(salmon_data_directory)
# create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE) {
dds <-
  create_dds('ZM', salmon_data_directory, times, "salmon_quant", 1:6, run_DESeq=FALSE)

annotation <- get_annotation()

# Where dds is a DESeqDataSet object
pcaExplorer(dds=dds,
            annotation=annotation)