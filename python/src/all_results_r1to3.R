


library("vsn")
library("genefilter")
library(dplyr)
library(openxlsx)

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
path <- "../../../Volumes/bs_external/lab_villunger"

 

load("r/data/ZM_r1to3_data.RData")
            
treatment = 'ZM'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/ZM/output_salmon'))
dds_ZM_r1to3 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_ZM_r1to3, file = glue('r/data/', "ZM_r1to3_data.RData"))

 

load("r/data/Noc_r1to3_data.RData")
            
treatment = 'Noc'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/Noc/output_salmon'))
dds_Noc_r1to3 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Noc_r1to3, file = glue('r/data/', "Noc_r1to3_data.RData"))

 

load("r/data/DHCB_r1to3_data.RData")
            
treatment = 'DHCB'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/DHCB/output_salmon'))
dds_DHCB_r1to3 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_DHCB_r1to3, file = glue('r/data/', "DHCB_r1to3_data.RData"))

 

load("r/data/Nutl_r1to3_data.RData")
            
treatment = 'Nutl'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/Nutl/output_salmon'))
dds_Nutl_r1to3 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Nutl_r1to3, file = glue('r/data/', "Nutl_r1to3_data.RData"))

 

load("r/data/Etop_r1to3_data.RData")
            
treatment = 'Etop'
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/Etop/output_salmon'))
dds_Etop_r1to3 <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_Etop_r1to3, file = glue('r/data/', "Etop_r1to3_data.RData"))



#### TIMEPOINT t0-t16

        
## TIMEPOINT t0-t16- loop for each time stamp
t0_t16_workbook <- createWorkbook()

# Create results 

addWorksheet(t0_t16_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_ZM_16 <- lfcShrink(dds_ZM_r1to3, coef="timepoint_t16_vs_t0", type="apeglm")
results_ZM_16 <- subset(results_ZM_16, padj < 0.1)  # Restrict to values which are significant
results_ZM_16 <- add_annotations_to_results(results_ZM_16)
results_ZM_16_df <- as.data.frame(results_ZM_16)

addWorksheet(t0_t16_workbook, "results_apelgm_ZM_16")
writeData(t0_t16_workbook, "results_subset_dataframe", results_ZM_16_df)

addWorksheet(t0_t16_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Noc_16 <- lfcShrink(dds_Noc_r1to3, coef="timepoint_t16_vs_t0", type="apeglm")
results_Noc_16 <- subset(results_Noc_16, padj < 0.1)  # Restrict to values which are significant
results_Noc_16 <- add_annotations_to_results(results_Noc_16)
results_Noc_16_df <- as.data.frame(results_Noc_16)

addWorksheet(t0_t16_workbook, "results_apelgm_Noc_16")
writeData(t0_t16_workbook, "results_subset_dataframe", results_Noc_16_df)

addWorksheet(t0_t16_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_DHCB_16 <- lfcShrink(dds_DHCB_r1to3, coef="timepoint_t16_vs_t0", type="apeglm")
results_DHCB_16 <- subset(results_DHCB_16, padj < 0.1)  # Restrict to values which are significant
results_DHCB_16 <- add_annotations_to_results(results_DHCB_16)
results_DHCB_16_df <- as.data.frame(results_DHCB_16)

addWorksheet(t0_t16_workbook, "results_apelgm_DHCB_16")
writeData(t0_t16_workbook, "results_subset_dataframe", results_DHCB_16_df)

addWorksheet(t0_t16_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Nutl_16 <- lfcShrink(dds_Nutl_r1to3, coef="timepoint_t16_vs_t0", type="apeglm")
results_Nutl_16 <- subset(results_Nutl_16, padj < 0.1)  # Restrict to values which are significant
results_Nutl_16 <- add_annotations_to_results(results_Nutl_16)
results_Nutl_16_df <- as.data.frame(results_Nutl_16)

addWorksheet(t0_t16_workbook, "results_apelgm_Nutl_16")
writeData(t0_t16_workbook, "results_subset_dataframe", results_Nutl_16_df)

addWorksheet(t0_t16_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Etop_16 <- lfcShrink(dds_Etop_r1to3, coef="timepoint_t16_vs_t0", type="apeglm")
results_Etop_16 <- subset(results_Etop_16, padj < 0.1)  # Restrict to values which are significant
results_Etop_16 <- add_annotations_to_results(results_Etop_16)
results_Etop_16_df <- as.data.frame(results_Etop_16)

addWorksheet(t0_t16_workbook, "results_apelgm_Etop_16")
writeData(t0_t16_workbook, "results_subset_dataframe", results_Etop_16_df)

saveWorkbook(t0_t16_workbook, "t0_t16_workbook.xlsx", overwrite = TRUE)



#### TIMEPOINT t0-t24

        
## TIMEPOINT t0-t24- loop for each time stamp
t0_t24_workbook <- createWorkbook()

# Create results 

addWorksheet(t0_t24_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_ZM_24 <- lfcShrink(dds_ZM_r1to3, coef="timepoint_t24_vs_t0", type="apeglm")
results_ZM_24 <- subset(results_ZM_24, padj < 0.1)  # Restrict to values which are significant
results_ZM_24 <- add_annotations_to_results(results_ZM_24)
results_ZM_24_df <- as.data.frame(results_ZM_24)

addWorksheet(t0_t24_workbook, "results_apelgm_ZM_24")
writeData(t0_t24_workbook, "results_subset_dataframe", results_ZM_24_df)

addWorksheet(t0_t24_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Noc_24 <- lfcShrink(dds_Noc_r1to3, coef="timepoint_t24_vs_t0", type="apeglm")
results_Noc_24 <- subset(results_Noc_24, padj < 0.1)  # Restrict to values which are significant
results_Noc_24 <- add_annotations_to_results(results_Noc_24)
results_Noc_24_df <- as.data.frame(results_Noc_24)

addWorksheet(t0_t24_workbook, "results_apelgm_Noc_24")
writeData(t0_t24_workbook, "results_subset_dataframe", results_Noc_24_df)

addWorksheet(t0_t24_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_DHCB_24 <- lfcShrink(dds_DHCB_r1to3, coef="timepoint_t24_vs_t0", type="apeglm")
results_DHCB_24 <- subset(results_DHCB_24, padj < 0.1)  # Restrict to values which are significant
results_DHCB_24 <- add_annotations_to_results(results_DHCB_24)
results_DHCB_24_df <- as.data.frame(results_DHCB_24)

addWorksheet(t0_t24_workbook, "results_apelgm_DHCB_24")
writeData(t0_t24_workbook, "results_subset_dataframe", results_DHCB_24_df)

addWorksheet(t0_t24_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Nutl_24 <- lfcShrink(dds_Nutl_r1to3, coef="timepoint_t24_vs_t0", type="apeglm")
results_Nutl_24 <- subset(results_Nutl_24, padj < 0.1)  # Restrict to values which are significant
results_Nutl_24 <- add_annotations_to_results(results_Nutl_24)
results_Nutl_24_df <- as.data.frame(results_Nutl_24)

addWorksheet(t0_t24_workbook, "results_apelgm_Nutl_24")
writeData(t0_t24_workbook, "results_subset_dataframe", results_Nutl_24_df)

addWorksheet(t0_t24_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Etop_24 <- lfcShrink(dds_Etop_r1to3, coef="timepoint_t24_vs_t0", type="apeglm")
results_Etop_24 <- subset(results_Etop_24, padj < 0.1)  # Restrict to values which are significant
results_Etop_24 <- add_annotations_to_results(results_Etop_24)
results_Etop_24_df <- as.data.frame(results_Etop_24)

addWorksheet(t0_t24_workbook, "results_apelgm_Etop_24")
writeData(t0_t24_workbook, "results_subset_dataframe", results_Etop_24_df)

saveWorkbook(t0_t24_workbook, "t0_t24_workbook.xlsx", overwrite = TRUE)



#### TIMEPOINT t0-t48

        
## TIMEPOINT t0-t48- loop for each time stamp
t0_t48_workbook <- createWorkbook()

# Create results 

addWorksheet(t0_t48_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_ZM_48 <- lfcShrink(dds_ZM_r1to3, coef="timepoint_t48_vs_t0", type="apeglm")
results_ZM_48 <- subset(results_ZM_48, padj < 0.1)  # Restrict to values which are significant
results_ZM_48 <- add_annotations_to_results(results_ZM_48)
results_ZM_48_df <- as.data.frame(results_ZM_48)

addWorksheet(t0_t48_workbook, "results_apelgm_ZM_48")
writeData(t0_t48_workbook, "results_subset_dataframe", results_ZM_48_df)

addWorksheet(t0_t48_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Noc_48 <- lfcShrink(dds_Noc_r1to3, coef="timepoint_t48_vs_t0", type="apeglm")
results_Noc_48 <- subset(results_Noc_48, padj < 0.1)  # Restrict to values which are significant
results_Noc_48 <- add_annotations_to_results(results_Noc_48)
results_Noc_48_df <- as.data.frame(results_Noc_48)

addWorksheet(t0_t48_workbook, "results_apelgm_Noc_48")
writeData(t0_t48_workbook, "results_subset_dataframe", results_Noc_48_df)

addWorksheet(t0_t48_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_DHCB_48 <- lfcShrink(dds_DHCB_r1to3, coef="timepoint_t48_vs_t0", type="apeglm")
results_DHCB_48 <- subset(results_DHCB_48, padj < 0.1)  # Restrict to values which are significant
results_DHCB_48 <- add_annotations_to_results(results_DHCB_48)
results_DHCB_48_df <- as.data.frame(results_DHCB_48)

addWorksheet(t0_t48_workbook, "results_apelgm_DHCB_48")
writeData(t0_t48_workbook, "results_subset_dataframe", results_DHCB_48_df)

addWorksheet(t0_t48_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Nutl_48 <- lfcShrink(dds_Nutl_r1to3, coef="timepoint_t48_vs_t0", type="apeglm")
results_Nutl_48 <- subset(results_Nutl_48, padj < 0.1)  # Restrict to values which are significant
results_Nutl_48 <- add_annotations_to_results(results_Nutl_48)
results_Nutl_48_df <- as.data.frame(results_Nutl_48)

addWorksheet(t0_t48_workbook, "results_apelgm_Nutl_48")
writeData(t0_t48_workbook, "results_subset_dataframe", results_Nutl_48_df)

addWorksheet(t0_t48_workbook, "t0_t16_workbook")

### loop for each treatment Template:
results_Etop_48 <- lfcShrink(dds_Etop_r1to3, coef="timepoint_t48_vs_t0", type="apeglm")
results_Etop_48 <- subset(results_Etop_48, padj < 0.1)  # Restrict to values which are significant
results_Etop_48 <- add_annotations_to_results(results_Etop_48)
results_Etop_48_df <- as.data.frame(results_Etop_48)

addWorksheet(t0_t48_workbook, "results_apelgm_Etop_48")
writeData(t0_t48_workbook, "results_subset_dataframe", results_Etop_48_df)

saveWorkbook(t0_t48_workbook, "t0_t48_workbook.xlsx", overwrite = TRUE)


