
library(DESeq2)
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

load('~/bs-villungerlab/results/output_encode_1to6/Noc_star_data.RData')
dds_Noc <- ddseq_Noc

results_Noc_t16_df <- return_results(dds_Noc, "timepoint_t16_vs_t0", "_16")
results_Noc_t20_df <- return_results(dds_Noc, "timepoint_t20_vs_t0", "_20")
results_Noc_t24_df <- return_results(dds_Noc, "timepoint_t24_vs_t0", "_24")
results_Noc_t36_df <- return_results(dds_Noc, "timepoint_t36_vs_t0", "_36")
results_Noc_t48_df <- return_results(dds_Noc, "timepoint_t48_vs_t0", "_48")
df <- results_Noc_t24_df

#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange), ]
# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_merged_df <- merge_all_data(upr_top, results_Noc_t16_df, results_Noc_t20_df, results_Noc_t36_df, results_Noc_t48_df, 'results/generegulation/Noc_gene_regulation_data.csv')
upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, 24)
upr_plot <- plot_longdf(upr_top_long_df, "Noc upregulated genes")
upr_plot
ggsave(filename = "results/generegulation/Noc_upregulated_genes.pdf", plot = upr_plot)

downr_df_sorted <- df[order(df$log2FoldChange), ]
# set the number of results which you want
downr_top <- head(upr_df_sorted, 25)
downr_top_merged_df <- merge_all_data(downr_top, dfs)
downr_top_long_df <- make_longdf_for_plot(downr_top_merged_df, 24)
downr_plot <- plot_longdf(upr_top_long_df, "Noc downregulated genes")
downr_plot
ggsave(filename = "results/generegulation/Noc_downregulated_genes.pdf", plot = downr_plot)
