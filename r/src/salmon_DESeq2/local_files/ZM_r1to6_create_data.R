
library("vsn")
library("openxlsx")
library("genefilter")
setwd("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

times = c(0, 16, 20, 24, 36, 48)
treatment <- "ZM"
data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_salmon'))


dds_ZM_r1to6 <- create_dds('ZM', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds_ZM_r1to6, file = glue('r/data/', "ZM_r1to6_data.RData"))

res_ZM_r1to6 <- results(dds_ZM_r1to6)
res_ZM_r1to6Ordered <- res_ZM_r1to6[order(res_ZM_r1to6$padj),]
res_ZM_r1to6Ordered <- add_annotations_to_results(res_ZM_r1to6Ordered)

head(res_ZM_r1to6Ordered)

res_ZM_r1to6OrderedDF <- as.data.frame(res_ZM_r1to6Ordered)
# write.csv(res_ZM_r1to6OrderedDF, file = "res_ZM_r1to6ults/ZM_r1to6_data.csv")

# Creating spreadsheets for analysis:::
zm_workbook <- createWorkbook()
addWorksheet(zm_workbook, "ZM_all")
writeData(zm_workbook, "ZM_all", res_ZM_r1to6OrderedDF)

res_vsd_ZM_r1to6 <- vst(dds_ZM_r1to6, blind=FALSE)
res_vsd_ZM_r1to6 <- results(dds_ZM_r1to6)
res_vsd_ZM_r1to6_ordered <- res_vsd_ZM_r1to6[order(res_vsd_ZM_r1to6$padj),]
res_vsd_ZM_r1to6_ordered <- add_annotations_to_results(res_vsd_ZM_r1to6_ordered)
res_vsd_ZM_r1to6_orderedDF <- as.data.frame(res_vsd_ZM_r1to6_ordered)

addWorksheet(zm_workbook, "ZM_vst")
writeData(zm_workbook, "ZM_vst", res_vsd_ZM_r1to6_orderedDF)

rld_ZM_r1to6 <- rlog(dds, blind=FALSE)
res_rld_ZM_r1to6 <- results(dds_ZM_r1to6)
res_rld_ZM_r1to6_ordered <- res_rld_ZM_r1to6[order(res_rld_ZM_r1to6$padj),]
res_rld_ZM_r1to6_ordered <- add_annotations_to_results(res_rld_ZM_r1to6_ordered)
res_rld_ZM_r1to6_orderedDF <- as.data.frame(res_rld_ZM_r1to6_ordered)
addWorksheet(zm_workbook, "ZM_rLog")
writeData(zm_workbook, "ZM_rLog", res_rld_ZM_r1to6_orderedDF)

shrinkage <- c("ashr", "apeglm")

# Save the workbook to an .xlsx file
saveWorkbook(zm_workbook, "zm_all_treatment_data.xlsx", overwrite = TRUE)

resultsNames(dds)
resAsh <- lfcShrink(dds_ZM_r1to6, coef='timepoint_t16_vs_t0', type="ashr")
