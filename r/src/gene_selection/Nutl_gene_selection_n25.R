#####
library(DESeq2)
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

setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")
######
load('~/bs-villungerlab/results/output_encode_1to6/Nutl_star_data.RData')

dpi = 500
width_in <- 12
height_in <- 12

dds_Nutl <- ddseq_Nutl
data_file = "~/bs-villungerlab/results/output_encode_1to6/Nutl_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_Nutl_t8_df <- return_results(dds_Nutl, "timepoint_t8_vs_t0", "_8")
  results_Nutl_t12_df <- return_results(dds_Nutl, "timepoint_t12_vs_t0", "_12")
  results_Nutl_t16_df <- return_results(dds_Nutl, "timepoint_t16_vs_t0", "_16")
  results_Nutl_t24_df <- return_results(dds_Nutl, "timepoint_t24_vs_t0", "_24")
  results_Nutl_t48_df <- return_results(dds_Nutl, "timepoint_t48_vs_t0", "_48")
  save(results_Nutl_t48_df, results_Nutl_t8_df, results_Nutl_t12_df, results_Nutl_t16_df, results_Nutl_t24_df, file=data_file)
}

df <- results_Nutl_t16_df

all_df_merged_df <- merge_all_data(results_Nutl_t48_df, results_Nutl_t8_df, results_Nutl_t12_df, results_Nutl_t16_df, results_Nutl_t24_df, 'results/output_encode/Nutl/all_Nutl_gene_regulation_data.csv', 'full_join')
df <- all_df_merged_df

df$symbol_48
df <- fix_labels(df)
df$symbol

#####
upr_df_sorted <- df[order(-df$log2FoldChange_16), ]
head(upr_df_sorted)
# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_long_df <- make_longdf_for_plot(upr_top, 16)
upr_plot <- plot_longdf(upr_top_long_df, "Nutl upregulated genes: n25 | t16")
upr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_16), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
downr_top_long_df <- make_longdf_for_plot(downr_top, 16)
downr_plot <- plot_longdf(downr_top_long_df, "Nutl downregulated genes: (3n) n25 | t16")
downr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)


#####
df 
dim(df)
# This join will conditionally affect the results, because it will restrict it based if the values are already there.
dim(all_df_merged_df)
all_df_merged_df$log2FoldChange_8
all_df_merged_df
#####

# subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_8))
middle_points <- subset(all_df_merged_df, !is.na(log2FoldChange_12))
middle_points <- subset(middle_points, !is.na(log2FoldChange_16))
middle_points <- subset(middle_points, !is.na(log2FoldChange_24))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_48))
dim(middle_points)
middle_points_long_df <- make_longdf_for_plot(middle_points, 16)

middle_points_plot <- plot_longdf(middle_points_long_df, "Nutl genes (12,16,24 timepoints)")
middle_points_plot

write.csv(middle_points, file = "results/output_encode/Nutl/Nutl_middle_timepoints_signature.csv")
ggsave(filename = "results/output_encode/Nutl/Nutl_middle_timepoints.pdf", plot = middle_points_plot, dpi=dpi, width=width_in, height=height_in)
######

full_signature <- subset(all_df_merged_df, !is.na(log2FoldChange_8))
full_signature <- subset(full_signature, !is.na(log2FoldChange_12))
full_signature <- subset(full_signature, !is.na(log2FoldChange_16))
full_signature <- subset(full_signature, !is.na(log2FoldChange_24))
full_signature <- subset(full_signature, !is.na(log2FoldChange_48))


full_signature$symbol_48
full_signature <- fix_labels(full_signature)
full_signature$symbol

                       
# write.csv(full_signature, file = "results/output_encode/Nutl/Nutl_5n_signature.csv")
full_signature_threshold <- full_signature[
  any(abs(full_signature$log2FoldChange_8)>1 | 
        abs(full_signature$log2FoldChange_12)>1 | 
        abs(full_signature$log2FoldChange_16)>1 | 
        abs(full_signature$log2FoldChange_24)>1 |
        abs(full_signature$log2FoldChange_48)>1), ]

write.csv(full_signature_threshold, file = "results/output_encode/Nutl/Nutl_5n_threshold_1lfc_signature.csv")

full_signature_threshold_decrease <- subset(full_signature_threshold, (log2FoldChange_16 < 0 | log2FoldChange_8 < 0 | log2FoldChange_12 < 0 | log2FoldChange_24 < 0 | log2FoldChange_48 < 0))
dim(full_signature_threshold_decrease)
write.csv(full_signature_threshold_decrease, file = "results/output_encode/Nutl/Nutl_full_signature_5n_decrease.csv")
noquote_full_signature_threshold_decrease_df <- noquote(full_signature_threshold_decrease$symbol)
write(noquote_full_signature_threshold_decrease_df, file = "results/output_encode/Nutl/Nutl_full_signature_5n_decrease.txt")

full_signature_threshold_increase <- subset(full_signature_threshold, (log2FoldChange_16 > 0 | log2FoldChange_8 > 0 | log2FoldChange_12 > 0 | log2FoldChange_24 > 0 | log2FoldChange_48 > 0))
dim(full_signature_threshold_increase)
write.csv(full_signature_threshold_increase, file = "results/output_encode/Nutl/Nutl_full_signature_5n_increase.csv")
noquote_full_signature_threshold_increase_df <- noquote(full_signature_threshold_increase$symbol)
write(noquote_full_signature_threshold_increase_df, file = "results/output_encode/Nutl/Nutl_full_signature_5n_increase.txt")

full_signature_threshold_long_df <- make_longdf_for_plot(full_signature_threshold, 16)
plot_title <- "Full Signature Nutl treatment Signature: l4c>1"
full_signature_plot <- ggplot(full_signature_threshold_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot")
full_signature_plot

ggsave(filename = "results/output_encode/Nutl/Nutl_5n_plot.pdf", plot = full_signature_plot, dpi=dpi, width=width_in, height=height_in)


df <- all_df_merged_df
df$symbol_48

df$symbol_48
df <- fix_labels(df)
df$symbol

# subset the data to only include data which is greater than 1 l2foldchange

subset_increase <- df[
  any(abs(df$log2FoldChange_8)>1 | 
  abs(df$log2FoldChange_12)>1 | 
  abs(df$log2FoldChange_16)>1 | 
  abs(df$log2FoldChange_24)>1 |
  abs(df$log2FoldChange_48)>1), ]

dim(df)
dim(subset_increase)

two_consecutivesubset_df <- subset_increase[
  !is.na(subset_increase$log2FoldChange_8) & !is.na(subset_increase$log2FoldChange_12) |
  !is.na(subset_increase$log2FoldChange_12) & !is.na(subset_increase$log2FoldChange_16) |
  !is.na(subset_increase$log2FoldChange_16) & !is.na(subset_increase$log2FoldChange_24) |
  !is.na(subset_increase$log2FoldChange_24) & !is.na(subset_increase$log2FoldChange_48)
  , ]

dim(df)
dim(two_consecutivesubset_df)

write.csv(two_consecutivesubset_df, file = "results/output_encode/Nutl/Nutl_gene_signature_2n.csv")
two_consecutivesubset_df_increase <- subset(two_consecutivesubset_df, (log2FoldChange_16 > 0 | log2FoldChange_8 > 0 | log2FoldChange_12 > 0 | log2FoldChange_24 > 0 | log2FoldChange_48 > 0))
dim(two_consecutivesubset_df_increase)
write.csv(two_consecutivesubset_df_increase, file = "results/output_encode/Nutl/Nutl_gene_signature_2n_increase.csv")

noquote_two_consecutivesubset_df_increase <- noquote(two_consecutivesubset_df_increase$symbol)
write(noquote_two_consecutivesubset_df_increase, file = "results/output_encode/Nutl/Nutl_gene_signature_2n_increase.txt")

two_consecutivesubset_df_decrease <- subset(two_consecutivesubset_df, (log2FoldChange_16 < 0 | log2FoldChange_8 < 0 | log2FoldChange_12 < 0 | log2FoldChange_24 < 0 | log2FoldChange_48 < 0))
dim(two_consecutivesubset_df_decrease)
write.csv(two_consecutivesubset_df_decrease, file = "results/output_encode/Nutl/Nutl_gene_signature_2n_decrease.csv")
noquote_two_consecutivesubset_df_decrease <- noquote(two_consecutivesubset_df_decrease$symbol)
write(noquote_two_consecutivesubset_df_decrease, file = "results/output_encode/Nutl/Nutl_gene_signature_2n_decrease.txt")

# write.table(df$column_name, file = "output.txt", row.names = FALSE, col.names = FALSE)

dim(two_consecutivesubset_df)
colnames(subset_increase)

three_consecutivesubset_df <- subset_increase[(!is.na(subset_increase$log2FoldChange_8) & !is.na(subset_increase$log2FoldChange_12) & !is.na(subset_increase$log2FoldChange_16) |
  !is.na(subset_increase$log2FoldChange_12) & !is.na(subset_increase$log2FoldChange_16) & !is.na(subset_increase$log2FoldChange_24) |
  !is.na(subset_increase$log2FoldChange_16) & !is.na(subset_increase$log2FoldChange_24) & !is.na(subset_increase$log2FoldChange_48)), ]

# dim(three_consecutivesubset_df)


write.csv(three_consecutivesubset_df, file = "results/output_encode/Nutl/Nutl_gene_signature_3n.csv")

three_consecutivesubset_df_increase <- subset(three_consecutivesubset_df, (log2FoldChange_16 > 0 | log2FoldChange_8 > 0 | log2FoldChange_12 > 0 | log2FoldChange_24 > 0 | log2FoldChange_48 > 0))
dim(three_consecutivesubset_df_increase)
write.csv(three_consecutivesubset_df_increase, file = "results/output_encode/Nutl/Nutl_gene_signature_3n_increase.csv")
noquote_three_consecutivesubset_df_increase <- noquote(three_consecutivesubset_df_increase$symbol)
write(noquote_three_consecutivesubset_df_increase, file = "results/output_encode/Nutl/Nutl_gene_signature_3n_increase.txt")

three_consecutivesubset_df_decrease <- subset(three_consecutivesubset_df, (log2FoldChange_16 < 0 | log2FoldChange_8 < 0 | log2FoldChange_12 < 0 | log2FoldChange_24 < 0 | log2FoldChange_48 < 0))
dim(three_consecutivesubset_df_decrease)
write.csv(three_consecutivesubset_df_decrease, file = "results/output_encode/Nutl/Nutl_gene_signature_3n_decrease.csv")
noquote_three_consecutivesubset_df_decrease <- noquote(three_consecutivesubset_df_decrease$symbol)
write(noquote_three_consecutivesubset_df_decrease, file = "results/output_encode/Nutl/Nutl_gene_signature_3n_decrease.txt")

dim(three_consecutivesubset_df)

three_consecutivesubset_df

two_consecutive_values_long_df <- make_longdf_for_plot(two_consecutivesubset_df, 16)
three_consecutivesubset_long_df <- make_longdf_for_plot(three_consecutivesubset_df, 16)
plot_title <- "2 Consecutive Values Nutl treatment Signature"
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

ggsave(filename = "results/output_encode/Nutl/Nutl_plot_two_consecutive_values.pdf", plot = two_p, dpi=dpi, width=width_in, height=height_in)


three_consecutivesubset_long_df <- make_longdf_for_plot(three_consecutivesubset_df, 16)
plot_title <- "3 Consecutive Values Nutl treatment Signature"
three_p <- ggplot(three_consecutivesubset_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "none")

ggsave(filename = "results/output_encode/Nutl/Nutl_plot_three_consecutive_values.pdf", plot = three_p, dpi=dpi, width=width_in, height=height_in)

three_p

df <- all_df_merged_df
#now you can apply filtering
df <- subset(df, !is.na(log2FoldChange_12))
df <- subset(df, !is.na(log2FoldChange_16))
df <- subset(df, !is.na(log2FoldChange_24))

upr_df_sorted <- df[order(-df$log2FoldChange_16), ]

# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_long_df <- make_longdf_for_plot(upr_top, 16)
upr_plot <- plot_longdf(upr_top_long_df, "Nutl upregulated genes: (3n) n25 | t16")
upr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_three_consecutive_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_16), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
downr_top_long_df <- make_longdf_for_plot(downr_top, 16)
downr_plot <- plot_longdf(downr_top_long_df, "Nutl downregulated genes: (3n) n25 | t16")
downr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_three_consecutive_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)

dim(df)
#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange_16), ]

# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)

df <- upr_top
df <- subset(df, !is.na(log2FoldChange_8))
df <- subset(df, !is.na(log2FoldChange_12))
df <- subset(df, !is.na(log2FoldChange_16))
df <- subset(df, !is.na(log2FoldChange_24))
df <- subset(df, !is.na(log2FoldChange_48))

upr_top_long_df <- make_longdf_for_plot(df, 16)
upr_plot <- plot_longdf(upr_top_long_df, "Nutl upregulated genes: (5n) n25 | t16")
upr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_full_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_16), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
df <- downr_top
df <- subset(df, !is.na(log2FoldChange_8))
df <- subset(df, !is.na(log2FoldChange_12))
df <- subset(df, !is.na(log2FoldChange_16))
df <- subset(df, !is.na(log2FoldChange_24))
df <- subset(df, !is.na(log2FoldChange_48))

downr_top_long_df <- make_longdf_for_plot(df, 16)
downr_plot <- plot_longdf(downr_top_long_df, "Nutl downregulated genes: (3n) n25 | t16")
downr_plot
ggsave(filename = "results/output_encode/Nutl/n25generegulation/Nutl_full_consecutive_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)
