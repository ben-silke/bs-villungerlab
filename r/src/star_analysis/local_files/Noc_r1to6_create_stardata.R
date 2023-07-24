
setwd("/Users/bsilke/bs-villungerlab")
library("glue")
library("Rsubread")
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
treatment <- "Noc"
data_directory = file.path('data/output_htseq_counts_2')
ddseq_Noc <- load_all_htseq_data(file.path(data_directory, 'all_Noc_fc.tsv'))

# <- create_htseq_ddseq(Noc, data_directory, times, 1:6)

save(ddseq_Noc, file = glue('r/data/', 'Noc_r1to6_star.RData'))

Noc_workbook <- createWorkbook()
times = c(16, 20, 24, 36, 48)
for (time in times) {
    timepoint <- glue("timepoint_t{time}_vs_t0")
    results_Noc <- lfcShrink(ddseq_Noc, coef=timepoint, type="apeglm")

    results_Noc <- subset(results_Noc, padj < 0.1)  # Restrict to values which are significant
    ## results_Noc <- add_annotations_to_results(results_Noc)
    results_Noc_df <- as.data.frame(results_Noc)
    addWorksheet(Noc_workbook, glue("Noc_{time}"))
    writeData(Noc_workbook, glue("Noc_{time}"), results_Noc_df)
}

saveWorkbook(Noc_workbook, "results/Noc_workbook.xlsx", overwrite = TRUE)

