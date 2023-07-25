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


setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_analysis/star_utils.R")
source("~/bs-villungerlab/r/src/star_analysis/gene_selection_utils.R")


load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')
dds_ZM <- ddseq_ZM

results_ZM_t24_df <- return_results(dds_ZM, "timepoint_t24_vs_t0", "_24")
results_ZM_t16_df <- return_results(dds_ZM, "timepoint_t16_vs_t0", "_16")
results_ZM_t20_df <- return_results(dds_ZM, "timepoint_t20_vs_t0", "_20")
results_ZM_t36_df <- return_results(dds_ZM, "timepoint_t36_vs_t0", "_36")
results_ZM_t48_df <- return_results(dds_ZM, "timepoint_t48_vs_t0", "_48")

df <- results_ZM_t24_df
# dfs
#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange),]
# set the number of results which you want
upr_top <- head(upr_df_sorted, 20)
upr_top
upr_top_merged_df <- merge_all_data(upr_top, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t36_df, results_ZM_t48_df, 'generegulation/ZM_gene_upregulation_data.csv' )
upr_top_merged_df
upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df)
upr_top_long_df
upr_plot <- plot_longdf(upr_top_long_df, "ZM upregulated genes")
upr_plot
ggsave(filename = "results/generegulation/ZM_upregulated_genes.pdf", plot = upr_plot)


downr_df_sorted <- df[order(df$log2FoldChange), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 20)
downr_top_merged_df <- merge_all_data(downr_top, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t36_df, results_ZM_t48_df, 'generegulation/ZM_gene_downregulation_data.csv')
downr_top_long_df <- make_longdf_for_plot(downr_top_merged_df)
downr_plot <- plot_longdf(downr_top_long_df, "ZM downregulated genes")
downr_plot
ggsave(filename = "results/generegulation/ZM_downregulated_genes.pdf", plot = downr_plot)

