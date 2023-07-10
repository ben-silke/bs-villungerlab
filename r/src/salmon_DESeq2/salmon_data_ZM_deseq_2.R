getwd()
source("r/src/salmon_DESeq2/utils.R")

treatment <- "ZM"
salmon_data_directory = file.path(getwd(), glue('data/organised/{treatment}/output_salmon'))
times = c(0, 16, 24, 36, 48)

print(salmon_data_directory)
# create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE) {
dds <-
  create_dds('ZM', salmon_data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "_data.RData")))


results <- results(dds)
resultsNames(dds)


results_t0_to_t16 <-
  results(dds, contrast = c("timepoint", "t16", "t0"))
results_t0_to_t16

results_t0_to_t24 <-
  results(
    dds,
    alpha = 0.05,
    lfcThreshold = 1,
    contrast = c("timepoint", "t24", "t0")
  )

results_t0_to_t24


summary(rest0t12_stringent)
plotMA(rest0t24)

