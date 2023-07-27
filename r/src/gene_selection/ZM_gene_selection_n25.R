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
source("~/bs-villungerlab/r/src/star_analysis/star_utils.R")
source("~/bs-villungerlab/r/src/star_analysis/gene_selection_utils.R")
######
load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')

dpi = 500
width_in <- 20
height_in <- 20

dds_ZM <- ddseq_ZM
data_file = "~/bs-villungerlab/results/output_encode_1to6/ZM_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_ZM_t16_df <- return_results(dds_ZM, "timepoint_t16_vs_t0", "_16")
  results_ZM_t20_df <- return_results(dds_ZM, "timepoint_t20_vs_t0", "_20")
  results_ZM_t24_df <- return_results(dds_ZM, "timepoint_t24_vs_t0", "_24")
  results_ZM_t36_df <- return_results(dds_ZM, "timepoint_t36_vs_t0", "_36")
  results_ZM_t48_df <- return_results(dds_ZM, "timepoint_t48_vs_t0", "_48")
  save(results_ZM_t48_df, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t24_df, results_ZM_t36_df, file=data_file)
}

df <- results_ZM_t24_df

all_df_merged_df <- merge_all_data(results_ZM_t48_df, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t24_df, results_ZM_t36_df, 'results/output_encode/ZM/all_ZM_gene_regulation_data.csv', 'full_join')
df <- all_df_merged_df
#####
upr_df_sorted <- df[order(-df$log2FoldChange_24), ]
head(upr_df_sorted)
# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_long_df <- make_longdf_for_plot(upr_top, 24)
upr_plot <- plot_longdf(upr_top_long_df, "ZM upregulated genes: n25 | t24")
upr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_24), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
downr_top_long_df <- make_longdf_for_plot(downr_top, 24)
downr_plot <- plot_longdf(downr_top_long_df, "ZM downregulated genes: (3n) n25 | t24")
downr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)


#####
df 
dim(df)
# This join will conditionally affect the results, because it will restrict it based if the values are already there.
dim(all_df_merged_df)
all_df_merged_df$log2FoldChange_16
all_df_merged_df
#####

# subset_df <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
middle_points <- subset(all_df_merged_df, !is.na(log2FoldChange_20))
middle_points <- subset(middle_points, !is.na(log2FoldChange_24))
middle_points <- subset(middle_points, !is.na(log2FoldChange_36))
# subset_df <- subset(subset_df, !is.na(log2FoldChange_48))
dim(middle_points)
middle_points_long_df <- make_longdf_for_plot(middle_points, 24)

middle_points_plot <- plot_longdf(middle_points_long_df, "ZM genes (20,24,36 timepoints)")
middle_points_plot

write.csv(middle_points, file = "results/output_encode/ZM/ZM_middle_timepoints_signature.csv")
ggsave(filename = "results/output_encode/ZM/ZM_middle_timepoints.pdf", plot = middle_points_plot, dpi=dpi, width=width_in, height=height_in)
######

full_signature <- subset(all_df_merged_df, !is.na(log2FoldChange_16))
full_signature <- subset(full_signature, !is.na(log2FoldChange_20))
full_signature <- subset(full_signature, !is.na(log2FoldChange_24))
full_signature <- subset(full_signature, !is.na(log2FoldChange_36))
full_signature <- subset(full_signature, !is.na(log2FoldChange_48))
                       
# write.csv(full_signature, file = "results/output_encode/ZM/ZM_5n_signature.csv")
full_signature_threshold <- full_signature[
  any(abs(full_signature$log2FoldChange_16)>1 | 
        abs(full_signature$log2FoldChange_20)>1 | 
        abs(full_signature$log2FoldChange_24)>1 | 
        abs(full_signature$log2FoldChange_36)>1 |
        abs(full_signature$log2FoldChange_48)>1), ]

write.csv(full_signature_threshold, file = "results/output_encode/ZM/ZM_5n_threshold_1lfc_signature.csv")

full_signature_threshold_decrease <- subset(full_signature_threshold, (log2FoldChange_24 < 0 | log2FoldChange_16 < 0 | log2FoldChange_20 < 0 | log2FoldChange_36 < 0 | log2FoldChange_48 < 0))
dim(full_signature_threshold_decrease)
write.csv(full_signature_threshold_decrease, file = "results/output_encode/ZM/ZM_full_signature_5n_decrease.csv")
noquote_full_signature_threshold_decrease_df <- noquote(full_signature_threshold_decrease$symbol_48)
write(noquote_full_signature_threshold_decrease_df, file = "results/output_encode/ZM/ZM_full_signature_5n_decrease.txt")

full_signature_threshold_increase <- subset(full_signature_threshold, (log2FoldChange_24 > 0 | log2FoldChange_16 > 0 | log2FoldChange_20 > 0 | log2FoldChange_36 > 0 | log2FoldChange_48 > 0))
dim(full_signature_threshold_increase)
write.csv(full_signature_threshold_increase, file = "results/output_encode/ZM/ZM_full_signature_5n_increase.csv")
noquote_full_signature_threshold_increase_df <- noquote(full_signature_threshold_increase$symbol_48)
write(noquote_full_signature_threshold_increase_df, file = "results/output_encode/ZM/ZM_full_signature_5n_increase.txt")

full_signature_threshold_long_df <- make_longdf_for_plot(full_signature_threshold, 24)
plot_title <- "Full Signature ZM treatment Signature: l4c>1"
full_signature_plot <- ggplot(full_signature_threshold_long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y = expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot")
full_signature_plot

ggsave(filename = "results/output_encode/ZM/ZM_5n_plot.pdf", plot = full_signature_plot, dpi=dpi, width=width_in, height=height_in)


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



# dim(df)
# dim(subset_increase)
# subset_increase <- df[abs(df$log2FoldChange_16)>1, ]
# subset_increase <- subset_increase[abs(subset_increase$log2FoldChange_20)>1, ]
# subset_increase <- subset_increase[abs(subset_increase$log2FoldChange_24)>1, ]
# subset_increase <- subset_increase[abs(subset_increase$log2FoldChange_36)>1, ]
# subset_increase <- subset_increase[abs(subset_increase$log2FoldChange_48)>1, ]

dim(df)
dim(subset_increase)

two_consecutivesubset_df <- subset_increase[
  !is.na(subset_increase$log2FoldChange_16) & !is.na(subset_increase$log2FoldChange_20) |
  !is.na(subset_increase$log2FoldChange_20) & !is.na(subset_increase$log2FoldChange_24) |
  !is.na(subset_increase$log2FoldChange_24) & !is.na(subset_increase$log2FoldChange_36) |
  !is.na(subset_increase$log2FoldChange_36) & !is.na(subset_increase$log2FoldChange_48)
  , ]

dim(df)
dim(two_consecutivesubset_df)

write.csv(two_consecutivesubset_df, file = "results/output_encode/ZM/ZM_gene_signature_2n.csv")
two_consecutivesubset_df_increase <- subset(two_consecutivesubset_df, (log2FoldChange_24 > 0 | log2FoldChange_16 > 0 | log2FoldChange_20 > 0 | log2FoldChange_36 > 0 | log2FoldChange_48 > 0))
dim(two_consecutivesubset_df_increase)
write.csv(two_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_2n_increase.csv")

noquote_two_consecutivesubset_df_increase <- noquote(two_consecutivesubset_df_increase$symbol_48)
write(noquote_two_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_2n_increase.txt")

two_consecutivesubset_df_decrease <- subset(two_consecutivesubset_df, (log2FoldChange_24 < 0 | log2FoldChange_16 < 0 | log2FoldChange_20 < 0 | log2FoldChange_36 < 0 | log2FoldChange_48 < 0))
dim(two_consecutivesubset_df_decrease)
write.csv(two_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_2n_decrease.csv")
noquote_two_consecutivesubset_df_decrease <- noquote(two_consecutivesubset_df_decrease$symbol_48)
write(noquote_two_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_2n_decrease.txt")

# write.table(df$column_name, file = "output.txt", row.names = FALSE, col.names = FALSE)

dim(two_consecutivesubset_df)
colnames(subset_increase)

three_consecutivesubset_df <- subset_increase[(!is.na(subset_increase$log2FoldChange_16) & !is.na(subset_increase$log2FoldChange_20) & !is.na(subset_increase$log2FoldChange_24) |
  !is.na(subset_increase$log2FoldChange_20) & !is.na(subset_increase$log2FoldChange_24) & !is.na(subset_increase$log2FoldChange_36) |
  !is.na(subset_increase$log2FoldChange_24) & !is.na(subset_increase$log2FoldChange_36) & !is.na(subset_increase$log2FoldChange_48)), ]

# dim(three_consecutivesubset_df)


write.csv(three_consecutivesubset_df, file = "results/output_encode/ZM/ZM_gene_signature_3n.csv")

three_consecutivesubset_df_increase <- subset(three_consecutivesubset_df, (log2FoldChange_24 > 0 | log2FoldChange_16 > 0 | log2FoldChange_20 > 0 | log2FoldChange_36 > 0 | log2FoldChange_48 > 0))
dim(three_consecutivesubset_df_increase)
write.csv(three_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_3n_increase.csv")
noquote_three_consecutivesubset_df_increase <- noquote(three_consecutivesubset_df_increase$symbol_48)
write(noquote_three_consecutivesubset_df_increase, file = "results/output_encode/ZM/ZM_gene_signature_3n_increase.txt")

three_consecutivesubset_df_decrease <- subset(three_consecutivesubset_df, (log2FoldChange_24 < 0 | log2FoldChange_16 < 0 | log2FoldChange_20 < 0 | log2FoldChange_36 < 0 | log2FoldChange_48 < 0))
dim(three_consecutivesubset_df_decrease)
write.csv(three_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_3n_decrease.csv")
noquote_three_consecutivesubset_df_decrease <- noquote(three_consecutivesubset_df_decrease$symbol_48)
write(noquote_three_consecutivesubset_df_decrease, file = "results/output_encode/ZM/ZM_gene_signature_3n_decrease.txt")

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

ggsave(filename = "results/output_encode/ZM/ZM_plot_two_consecutive_values.pdf", plot = two_p, dpi=dpi, width=width_in, height=height_in)


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

ggsave(filename = "results/output_encode/ZM/ZM_plot_three_consecutive_values.pdf", plot = three_p, dpi=dpi, width=width_in, height=height_in)

three_p

df <- all_df_merged_df
#now you can apply filtering
df <- subset(df, !is.na(log2FoldChange_20))
df <- subset(df, !is.na(log2FoldChange_24))
df <- subset(df, !is.na(log2FoldChange_36))

upr_df_sorted <- df[order(-df$log2FoldChange_24), ]

# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)
upr_top_long_df <- make_longdf_for_plot(upr_top, 24)
upr_plot <- plot_longdf(upr_top_long_df, "ZM upregulated genes: (3n) n25 | t24")
upr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_three_consecutive_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_24), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
downr_top_long_df <- make_longdf_for_plot(downr_top, 24)
downr_plot <- plot_longdf(downr_top_long_df, "ZM downregulated genes: (3n) n25 | t24")
downr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_three_consecutive_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)

dim(df)
#now you can apply filtering
upr_df_sorted <- df[order(-df$log2FoldChange_24), ]

# set the number of results which you want
upr_top <- head(upr_df_sorted, 25)

df <- upr_top
df <- subset(df, !is.na(log2FoldChange_16))
df <- subset(df, !is.na(log2FoldChange_20))
df <- subset(df, !is.na(log2FoldChange_24))
df <- subset(df, !is.na(log2FoldChange_36))
df <- subset(df, !is.na(log2FoldChange_48))

upr_top_long_df <- make_longdf_for_plot(df, 24)
upr_plot <- plot_longdf(upr_top_long_df, "ZM upregulated genes: (5n) n25 | t24")
upr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_full_upregulated_genes.pdf", plot = upr_plot, dpi=dpi, width=width_in, height=height_in)

downr_df_sorted <- df[order(df$log2FoldChange_24), ]
# set the number of results which you want
downr_top <- head(downr_df_sorted, 25)
df <- downr_top
df <- subset(df, !is.na(log2FoldChange_16))
df <- subset(df, !is.na(log2FoldChange_20))
df <- subset(df, !is.na(log2FoldChange_24))
df <- subset(df, !is.na(log2FoldChange_36))
df <- subset(df, !is.na(log2FoldChange_48))

downr_top_long_df <- make_longdf_for_plot(df, 24)
downr_plot <- plot_longdf(downr_top_long_df, "ZM downregulated genes: (3n) n25 | t24")
downr_plot
ggsave(filename = "results/output_encode/ZM/n25generegulation/ZM_full_consecutive_downregulated_genes.pdf", plot = downr_plot, dpi=dpi, width=width_in, height=height_in)
