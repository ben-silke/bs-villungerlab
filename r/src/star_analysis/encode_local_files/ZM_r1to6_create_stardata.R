
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

source("r/src/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

times = c(0, 16, 20, 24, 36, 48)
treatment <- "ZM"
data_directory = file.path('data/encode_output_htseq_counts')
ddseq_ZM <- load_all_htseq_data(file.path(data_directory, 'all_ZM_htseq_encode_counts.tsv'))

# <- create_htseq_ddseq(ZM, data_directory, times, 1:6)

save(ddseq_ZM, file = 'results/output_encode_1to6/ZM_star_data.Rdata')

ZM_workbook <- createWorkbook()
times = c(16, 20, 24, 36, 48)
for (time in times) {
    timepoint <- glue("timepoint_t{time}_vs_t0")
    results_ZM <- lfcShrink(ddseq_ZM, coef=timepoint, type="apeglm")

    results_ZM <- subset(results_ZM, padj < 0.1)  # Restrict to values which are significant
    results_ZM <- add_annotations_to_results(results_ZM)
    results_ZM_df <- as.data.frame(results_ZM)
    addWorksheet(ZM_workbook, glue("ZM_{time}"))
    writeData(ZM_workbook, glue("ZM_{time}"), results_ZM_df, rowNames=TRUE)
}

saveWorkbook(ZM_workbook, "results/output_encode_1to6/ZM_workbook.xlsx", overwrite = TRUE)

