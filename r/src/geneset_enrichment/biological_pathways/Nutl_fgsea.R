#####
library(DESeq2)
library("apeglm")
library(org.Hs.eg.db)
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(htmlwidgets)
library("fgsea")
library(msigdbr)


setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")
######
load('~/bs-villungerlab/results/output_encode_1to6/Nutl_star_data.RData')
dir.create("results/output_encode/Nutl/fgsea_enrichment")


data_file = "~/bs-villungerlab/results/output_encode_1to6/Nutl_unfiltered_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_Nutl_t8 <- get_unfiltered_results(dds_Nutl, "timepoint_t8_vs_t0", "_8")
  results_Nutl_t12 <- get_unfiltered_results(dds_Nutl, "timepoint_t12_vs_t0", "_12")
  results_Nutl_t16 <- get_unfiltered_results(dds_Nutl, "timepoint_t16_vs_t0", "_16")
  results_Nutl_t24 <- get_unfiltered_results(dds_Nutl, "timepoint_t24_vs_t0", "_24")
  results_Nutl_t48 <- get_unfiltered_results(dds_Nutl, "timepoint_t48_vs_t0", "_48")
  save(results_Nutl_t48, results_Nutl_t8, results_Nutl_t12, results_Nutl_t16, results_Nutl_t24, file=data_file)
}

all_df_merged_df <- merge_all_data(results_Nutl_t48, results_Nutl_t8, results_Nutl_t12, results_Nutl_t16, results_Nutl_t24, 'results/output_encode/Nutl/all_Nutl_gene_regulation_data.csv', 'full_join')
df <- fix_labels(all_df_merged_df)

data_file = "~/bs-villungerlab/results/output_encode_1to6/misgdrbr_df.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  msigdbr_df = msigdbr(species = "human")
  save(msigdbr_df, file=data_file)
}

msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

abs_foldchange_increase <- df[
  any(abs(df$log2FoldChange_8)>1 | 
  abs(df$log2FoldChange_12)>1 | 
  abs(df$log2FoldChange_16)>1 | 
  abs(df$log2FoldChange_24)>1 |
  abs(df$log2FoldChange_48)>1), ]


# include only 3 consecutive timepoints.
subset_df_3n <- abs_foldchange_increase[(!is.na(abs_foldchange_increase$log2FoldChange_8) & !is.na(abs_foldchange_increase$log2FoldChange_12) & !is.na(abs_foldchange_increase$log2FoldChange_16) |
  !is.na(abs_foldchange_increase$log2FoldChange_12) & !is.na(abs_foldchange_increase$log2FoldChange_16) & !is.na(abs_foldchange_increase$log2FoldChange_24) |
  !is.na(abs_foldchange_increase$log2FoldChange_16) & !is.na(abs_foldchange_increase$log2FoldChange_24) & !is.na(abs_foldchange_increase$log2FoldChange_48)), ]

# Create Upper
######
# Need to sort the dataset before we provide to fgsea otherwise it will break.
sorted_subset_foldchange_increase_3n <- subset_df_3n[order(-subset_df_3n$log2FoldChange_16), ]

sorted_foldchange_increase_3n_tovec <- data.frame(gene_symbol = sorted_subset_foldchange_increase_3n$symbol, log2foldchange <- sorted_subset_foldchange_increase_3n$log2FoldChange_16)
sorted_foldchange_increase_3n_tovec <- subset(sorted_foldchange_increase_3n_tovec, !is.na(log2foldchange))
nvec_increase_3n <- setNames(sorted_foldchange_increase_3n_tovec$log2foldchange, sorted_foldchange_increase_3n_tovec$gene_symbol)

fgseaRes_increase_3n <- fgsea(pathways = msigdbr_list,
                  stats    = nvec_increase_3n,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes_increase_3n$leadingEdge <- sapply(fgseaRes_increase_3n$leadingEdge, toString)
write.csv(fgseaRes_increase_3n, file = "results/output_encode/Nutl/fgsea_enrichment/fgsea_enrichment_Nutl_increase_3n.csv")