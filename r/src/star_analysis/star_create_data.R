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

source("r/src/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")



times = c(0, 16, 20, 24, 36, 48)
treatment <- "ZM"
data_directory = file.path('/Volumes/bs_external/lab_villunger/htseq_counts/output_htseq_counts')

ddseq_ZM <- create_htseq_ddseq(treatment, data_directory, times, 1:6)

save(ddseq_ZM, file = glue('r/data/', "ZM_r1to6_STAR.RData"))

ZM_workbook <- createWorkbook() 

times = c(16, 20, 24, 36, 48)
ddseq_ZM

for (time in times) {
    results_ZM <- lfcShrink(ddseq_ZM, coef="timepoint_t16_vs_t0", type="apeglm")

    results_ZM <- subset(results_ZM, padj < 0.1)  # Restrict to values which are significant
    # results_ZM <- add_annotations_to_results(results_ZM)
    results_ZM_df <- as.data.frame(results_ZM)
    addWorksheet(ZM_workbook, glue("ZM_{time}"))
    writeData(ZM_workbook, glue("ZM_{time}"), results_ZM_df)
}
