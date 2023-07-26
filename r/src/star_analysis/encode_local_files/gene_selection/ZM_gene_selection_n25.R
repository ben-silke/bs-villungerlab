#####
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
######
load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')
dds_ZM <- ddseq_ZM

results_ZM_t16_df <- return_results(dds_ZM, "timepoint_t16_vs_t0", "_16")
results_ZM_t20_df <- return_results(dds_ZM, "timepoint_t20_vs_t0", "_20")
results_ZM_t24_df <- return_results(dds_ZM, "timepoint_t24_vs_t0", "_24")
results_ZM_t36_df <- return_results(dds_ZM, "timepoint_t36_vs_t0", "_36")
results_ZM_t48_df <- return_results(dds_ZM, "timepoint_t48_vs_t0", "_48")
df <- results_ZM_t24_df

#####
#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange), ]

# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_merged_df <- merge_all_data(upr_top, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t36_df, results_ZM_t48_df, 'results/output_encode/ZM/n25generegulation/ZM_gene_regulation_data.csv')
upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, 24)
upr_plot <- plot_longdf(upr_top_long_df, "ZM upregulated genes: n25 | t24")
upr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_upregulated_genes.pdf", plot = upr_plot)

downr_df_sorted <- df[order(df$log2FoldChange), ]
# set the number of results which you want
downr_top <- head(upr_df_sorted, 25)
downr_top_merged_df <- merge_all_data(downr_top, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t36_df, results_ZM_t48_df, 'results/output_encode/ZM/n25generegulation/ZM_gene_regulation_data.csv')
downr_top_long_df <- make_longdf_for_plot(downr_top_merged_df, 24)
downr_plot <- plot_longdf(downr_top_long_df, "ZM downregulated genes: n25 | t24")
downr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_downregulated_genes.pdf", plot = downr_plot)

#####
df 
dim(df)
# This join will conditionally affect the results, because it will restrict it based if the values are already there.
all_df_merged_df <- merge_all_data(results_ZM_t48_df, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t24_df, results_ZM_t36_df, 'results/output_encode/ZM/all_ZM_gene_regulation_data.csv', 'full_join')
dim(all_df_merged_df)
all_df_merged_df$log2FoldChange_16


#######
# subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_20))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_36))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_48))
# 
# # 'results/output_encode/ZM/all_ZM_gene_regulation_data.csv'
# write.csv(subset_df, file = "results/output_encode/ZM/ZM_all_timepoints.csv")
# subset_long_df <- make_longdf_for_plot(subset_df, 24)
# subset_plot <- plot_longdf(subset_long_df, "ZM genes (all timepoints)")
# subset_plot

# ggsave(filename = "results/output_encode/ZM/ZM_all_timepoints.pdf", plot = subset_plot)

# upr_top_long_df <- make_longdf_for_plot(upr_top_merged_df, 24)


all_df_merged_df
#####

# subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
semi_subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_20))
semi_subset_df <- subset(semi_subset_df, !is.na(log2FoldChange_24))
semi_subset_df <- subset(semi_subset_df, !is.na(log2FoldChange_36))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_48))
dim(semi_subset_df)
semisubset_long_df <- make_longdf_for_plot(semi_subset_df, 24)

semi_subset_plot <- plot_longdf(semisubset_long_df, "ZM genes (20,24,36 timepoints)")
semi_subset_plot

write.csv(semi_subset_df, file = "results/output_encode/ZM/ZM_middle_timepoints_signature.csv")
ggsave(filename = "results/output_encode/ZM/ZM_middle_timepoints.pdf", plot = semi_subset_plot)
######

full_signature <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
full_signature <- subset(full_signature, !is.na(log2FoldChange_20))
full_signature <- subset(full_signature, !is.na(log2FoldChange_24))
full_signature <- subset(full_signature, !is.na(log2FoldChange_36))
full_signature <- subset(full_signature, !is.na(log2FoldChange_48))
                       
write.csv(full_signature, file = "results/output_encode/ZM/ZM_full_signature_signature.csv")


df <- all_df_merged_df
df
dim(df)
colnames(df)
# subset the data to only include data which is greater than 1 l2foldchange
subset_increase <- df[
  any(abs(df$log2FoldChange_16)>1 | 
  abs(df$log2FoldChange_20)>1 | 
  abs(df$log2FoldChange_24)>1 | 
  abs(df$log2FoldChange_36)>1 |
  abs(df$log2FoldChange_48)>1), ]

dim(df)
dim(subset_increase)
# subset_df <- df[apply(df, 1, function(x) any(x["col1"] > 5, x["col2"] < 15)), ]

df <- subset_increase

two_consecutivesubset_df <- df[
  !is.na(df$log2FoldChange_16) & !is.na(df$log2FoldChange_20) |
  !is.na(df$log2FoldChange_20) & !is.na(df$log2FoldChange_24) |
  !is.na(df$log2FoldChange_24) & !is.na(df$log2FoldChange_36) |
  !is.na(df$log2FoldChange_36) & !is.na(df$log2FoldChange_48)
  , ]

write.csv(two_consecutivesubset_df, file = "results/output_encode/ZM/ZM_gene_signature_two_points.csv")
two_consecutivesubset_df_increase <- subset(two_consecutivesubset_df, (log2FoldChange_24 > 0 | log2FoldChange_16 > 0 | log2FoldChange_20 > 0 | log2FoldChange_36 > 0 | log2FoldChange_48 > 0))
dim(two_consecutivesubset_df_increase)
write.csv(two_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_two_points_increase.csv")

two_consecutivesubset_df_decrease <- subset(two_consecutivesubset_df, (log2FoldChange_24 < 0 | log2FoldChange_16 < 0 | log2FoldChange_20 < 0 | log2FoldChange_36 < 0 | log2FoldChange_48 < 0))
dim(two_consecutivesubset_df_decrease)
write.csv(two_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_two_points_decrease.csv")

dim(two_consecutivesubset_df)

three_consecutivesubset_df <- df[
  !is.na(df$log2FoldChange_16) & !is.na(df$log2FoldChange_20) & !is.na(df$log2FoldChange_24) |
    !is.na(df$log2FoldChange_20) & !is.na(df$log2FoldChange_24) & !is.na(df$log2FoldChange_36) |
    !is.na(df$log2FoldChange_24) & !is.na(df$log2FoldChange_36) & !is.na(df$log2FoldChange_48)
  , ]

write.csv(three_consecutivesubset_df, file = "results/output_encode/ZM/ZM_gene_signature_three_points.csv")

three_consecutivesubset_df_increase <- subset(three_consecutivesubset_df, (log2FoldChange_24 > 0 | log2FoldChange_16 > 0 | log2FoldChange_20 > 0 | log2FoldChange_36 > 0 | log2FoldChange_48 > 0))
dim(three_consecutivesubset_df_increase)
write.csv(three_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_three_points_increase.csv")

three_consecutivesubset_df_increase <- subset(three_consecutivesubset_df, (log2FoldChange_24 < 0 | log2FoldChange_16 < 0 | log2FoldChange_20 < 0 | log2FoldChange_36 < 0 | log2FoldChange_48 < 0))
dim(three_consecutivesubset_df_decrease)
write.csv(three_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_three_points_decrease.csv")

dim(three_consecutivesubset_df)
three_consecutivesubset_df

two_consecutive_values_long_df <- make_longdf_for_plot(two_consecutivesubset_df, 24)
three_consecutivesubset_long_df <- make_longdf_for_plot(three_consecutivesubset_df, 24)
plot_title <- "2 Consecutive Values ZM treatment Signature"
two_p <- ggplot(two_consecutive_values_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "none")
two_p

ggsave(filename = "results/output_encode/ZM/ZM_plot_two_consecutive_values.pdf", plot = two_p)


three_consecutivesubset_long_df <- make_longdf_for_plot(three_consecutivesubset_df, 24)
plot_title <- "3 Consecutive Values ZM treatment Signature"
three_p <- ggplot(three_consecutivesubset_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "none")

ggsave(filename = "results/output_encode/ZM/ZM_plot_three_consecutive_values.pdf", plot = three_p)

three_p
