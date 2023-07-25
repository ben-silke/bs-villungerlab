library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(tidyverse)

source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")
source("r/src/star_analysis/star_utils.R")

load('~/bs-villungerlab/results/output_encode_1to6/Noc_star_data.RData')
dds_Noc <- ddseq_Noc

results_Noc_t24_df <- return_results(dds, "timepoint_t24_vs_t0", "_24")
results_Noc_t16_df <- return_results(dds, "timepoint_t16_vs_t0", "_16")
results_Noc_t20_df <- return_results(dds, "timepoint_t20_vs_t0", "_20")
results_Noc_t36_df <- return_results(dds, "timepoint_t36_vs_t0", "_36")
results_Noc_t48_df <- return_results(dds, "timepoint_t48_vs_t0", "_48")

dfs <- list(results_Noc_t16_df, results_Noc_t20_df, results_Noc_t36_df, results_Noc_t48_df)
df <- results_Noc_t24_df

#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange), ]
upr_top <- head(upr_df_sorted, 20)
upr_top_merged_df <- merge_all_data(upr_top, dfs)
upr_top_long_df <- merge_all_data(upr_top_merged_df, dfs)
upr_plot <- plot_longdf(upr_top_long_df, "Noc upregulated genes")
upr_plot


downr_df_sorted <- df[order(df$log2FoldChange), ]
downr_top <- head(upr_df_sorted, 20)
downr_top_merged_df <- merge_all_data(downr_top, dfs)
downr_top_long_df <- merge_all_data(downr_top_merged_df, dfs)
downr_plot <- plot_longdf(upr_top_long_df, "Noc downregulated genes")
downr_plot

