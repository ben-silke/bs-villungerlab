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

getwd()
treatment = 'DHCB'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))
dds_DHCB_r1to6 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_DHCB_r1to6, file = glue('r/data/', "DHCB_r1to6_data.RData"))

treatment = 'Etop'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))
dds_Etop_r1to6 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Etop_r1to6, file = glue('r/data/', "Etop_r1to6_data.RData"))

treatment = 'Noc'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))
dds_Noc_r1to6 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Noc_r1to6, file = glue('r/data/', "Noc_r1to6_data.RData"))

treatment = 'Nutl'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))
dds_Nutl_r1to6 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Nutl_r1to6, file = glue('r/data/', "Nutl_r1to6_data.RData"))

treatment = 'ZM'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))
dds_ZM_r1to6 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_ZM_r1to6, file = glue('r/data/', "ZM_r1to6_data.RData"))



## TO COMPARE ALL:
#### TIMEPOINT t0-t16- loop for each time stamp
t0_t16_workbook <- createWorkbook()
addWorksheet(t0_t16_workbook, "t0_t16_workbook")

# Create results 
### loop for each treatment Template:
results <- lfcShrink(dds, coef=timepoint, type="apeglm")
results_subset <- subset(results, padj < 0.1)  # Restrict to values which are significant
results_subset <- add_annotations_to_results(results_subset)
results_subset_dataframe <- as.data.frame(results_subset)

addWorksheet(t0_t16_workbook, "results_subset_dataframe")
writeData(t0_t16_workbook, "results_subset_dataframe", results_subset_dataframe)


# Save worksheet at the end
saveWorkbook(zm_workbook, "zm_all_treatment_data.xlsx", overwrite = TRUE)

#### TIMEPOINT t0-t24


#### TIMEPOINT t0-t48



