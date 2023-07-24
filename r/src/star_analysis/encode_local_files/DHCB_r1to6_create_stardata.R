
setwd("/Users/bsilke/bs-villungerlab")
library("glue")
library("stringr")
library("DESeq2")
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library(dplyr)
library(openxlsx)

source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

times = c(0, 16, 20, 24, 36, 48)
treatment <- "DHCB"
data_directory = file.path('data/encode_output_htseq_counts')
ddseq_DHCB <- load_all_htseq_data(file.path(data_directory, 'all_DHCB_fc.tsv'))

# <- create_htseq_ddseq(DHCB, data_directory, times, 1:6)

save(ddseq_DHCB, file = 'results/output_encode_1to6/ZM_r1to6_star.RData')

DHCB_workbook <- createWorkbook()
times = c(16, 20, 24, 36, 48)
for (time in times) {
    timepoint <- glue("timepoint_t{time}_vs_t0")
    results_DHCB <- lfcShrink(ddseq_DHCB, coef=timepoint, type="apeglm")

    results_DHCB <- subset(results_DHCB, padj < 0.1)  # Restrict to values which are significant
    results_DHCB <- add_annotations_to_results(results_DHCB)
    results_DHCB_df <- as.data.frame(results_DHCB)
    addWorksheet(DHCB_workbook, glue("DHCB_{time}"))
    writeData(DHCB_workbook, glue("DHCB_{time}"), results_DHCB_df, row.names=TRUE)
}

saveWorkbook(DHCB_workbook, "results/output_encode_1to6/DHCB_workbook.xlsx", overwrite = TRUE)

