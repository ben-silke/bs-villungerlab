
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

# We can update whatever noc is
df <- results_Noc_t24_df



#######
#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange), ]
# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_merged_df <- merge_all_data(upr_top, results_Noc_t16_df, results_Noc_t20_df, results_Noc_t36_df, results_Noc_t48_df, 'results/output_encode/n25generegulation/Noc_gene_regulation_data.csv')
upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, 24)
upr_plot <- plot_longdf(upr_top_long_df, "Noc upregulated genes: n25 | t24")
upr_plot
ggsave(filename = "results/output_encode/n25generegulation/Noc_upregulated_genes.pdf", plot = upr_plot)

downr_df_sorted <- df[order(df$log2FoldChange), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
downr_top_merged_df <- merge_all_data(downr_top, results_Noc_t16_df, results_Noc_t20_df, results_Noc_t36_df, results_Noc_t48_df, 'results/output_encode/n25generegulation/Noc_gene_regulation_data.csv')
downr_top_long_df <- make_longdf_for_plot(downr_top_merged_df, 24)
downr_plot <- plot_longdf(downr_top_long_df, "Noc downregulated genes: n25 | t24")
downr_plot
ggsave(filename = "results/output_encode/n25generegulation/Noc_downregulated_genes.pdf", plot = downr_plot)

df 
all_df_merged_df <- merge_all_data(df, results_Noc_t16_df, results_Noc_t20_df, results_Noc_t36_df, results_Noc_t48_df, 'results/output_encode/all_Noc_gene_regulation_data.csv')
all_df_merged_df$log2FoldChange_16

subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
subset_df <- subset(subset_df, !is.na(log2FoldChange_20))
subset_df <- subset(subset_df, !is.na(log2FoldChange_36))
subset_df <- subset(subset_df, !is.na(log2FoldChange_48))

# 'results/output_encode/all_Noc_gene_regulation_data.csv'
write.csv(subset_df, file = "results/output_encode/Noc_all_timepoints.csv")

subset_long_df <- make_longdf_for_plot(subset_df, 24)
subset_plot <- plot_longdf(subset_long_df, "Noc genes (all timepoints)")
subset_plot

ggsave(filename = "results/output_encode/Noc_all_timepoints.pdf", plot = subset_plot)

# upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, 24)

all_df_merged_df

# subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
semi_subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_20))
semi_subset_df <- subset(semi_subset_df, !is.na(log2FoldChange_36))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_48))
dim(semi_subset_df)
semisubset_long_df <- make_longdf_for_plot(semi_subset_df, 24)

semi_subset_plot <- plot_longdf(semisubset_long_df, "Noc genes (20,24,36 timepoints)")
semi_subset_plot

write.csv(semi_subset_df, file = "results/output_encode/Noc_middle_timepoints.csv")
ggsave(filename = "results/output_encode/Noc_middle_timepoints.pdf", plot = semi_subset_plot)


df <- all_df_merged_df
df
dim(df)
colnames(df)
# subset the data to only include data which is greater than 1 l2foldchange
subset_increase <- df[
  abs(df$log2FoldChange_16)>1 | 
    abs(df$log2FoldChange_20)>1 | 
    abs(df$log2FoldChange_24)>1 | 
    abs(df$log2FoldChange_36)>1 |
    abs(df$log2FoldChange_48)>1, ]

subset_increase

# subset_df <- df[apply(df, 1, function(x) any(x["col1"] > 5, x["col2"] < 15)), ]

# Subset if two values are consecutive.
subset_df <- df[
  !is.na(df$log2FoldChange_16) & !is.na(df$log2FoldChange_20) |
    !is.na(df$log2FoldChange_20) & !is.na(df$log2FoldChange_24) |
    !is.na(df$log2FoldChange_24) & !is.na(df$log2FoldChange_36) |
    !is.na(df$log2FoldChange_36) & !is.na(df$log2FoldChange_48)
  , ]

two_consecutive_values_long_df <- make_longdf_for_plot(subset_df, 24)
plot_title <- "Two Consecutive Values Noc treatment Signature"
p <- ggplot(two_consecutive_values_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "bottom")
p
subset_increase


