
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
treatment <- "Etop"
data_directory = file.path('data/output_htseq_counts_2')
ddseq_Etop <- load_all_htseq_data(file.path(data_directory, 'all_Etop_fc.tsv'))

# <- create_htseq_ddseq(Etop, data_directory, times, 1:6)

save(ddseq_Etop, file = 'r/data/Etop_r1to6_star.RData')

Etop_workbook <- createWorkbook()
times = c(8, 12, 16, 24, 48)
for (time in times) {
    timepoint <- glue("timepoint_t{time}_vs_t0")
    results_Etop <- lfcShrink(ddseq_Etop, coef=timepoint, type="apeglm")

    results_Etop <- subset(results_Etop, padj < 0.1)  # Restrict to values which are significant
    ## results_Etop <- add_annotations_to_results(results_Etop)
    results_Etop_df <- as.data.frame(results_Etop)
    addWorksheet(Etop_workbook, glue("Etop_{time}"))
    writeData(Etop_workbook, glue("Etop_{time}"), results_Etop_df)
}

saveWorkbook(Etop_workbook, "results/Etop_workbook.xlsx", overwrite = TRUE)

