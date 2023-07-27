
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

times = c(0, 8, 12, 16, 24, 48)
treatment <- "Nutl"
data_directory = file.path('data/output_htseq_counts_2')
ddseq_Nutl <- load_all_htseq_data(file.path(data_directory, 'all_Nutl_fc.tsv'))

# <- create_htseq_ddseq(Nutl, data_directory, times, 1:6)

save(ddseq_Nutl, file = 'r/data/Nutl_r1to6_star.RData')

Nutl_workbook <- createWorkbook()
times = c(8, 12, 16, 24, 48)
for (time in times) {
    timepoint <- glue("timepoint_t{time}_vs_t0")
    results_Nutl <- lfcShrink(ddseq_Nutl, coef=timepoint, type="apeglm")

    results_Nutl <- subset(results_Nutl, padj < 0.1)  # Restrict to values which are significant
    results_Nutl <- add_annotations_to_results(results_Nutl)
    results_Nutl_df <- as.data.frame(results_Nutl)
    addWorksheet(Nutl_workbook, glue("Nutl_{time}"))
    writeData(Nutl_workbook, glue("Nutl_{time}"), results_Nutl_df)
}

saveWorkbook(Nutl_workbook, "results/Nutl_workbook.xlsx", overwrite = TRUE)

