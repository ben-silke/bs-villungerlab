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


load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')

dds_zm <- ddseq_ZM

results_zm_t24 <- lfcShrink(dds_zm, coef="timepoint_t24_vs_t0", type="apeglm")
results_zm_t24Sig <- subset(results_zm_t24, padj < 0.1)
results_zm_t24Sig <- add_annotations_to_results(results_zm_t24Sig)

results_zm_t16 <- lfcShrink(dds_zm, coef="timepoint_t16_vs_t0", type="apeglm")
results_zm_t20 <- lfcShrink(dds_zm, coef="timepoint_t20_vs_t0", type="apeglm")
results_zm_t36 <- lfcShrink(dds_zm, coef="timepoint_t36_vs_t0", type="apeglm")
results_zm_t48 <- lfcShrink(dds_zm, coef="timepoint_t48_vs_t0", type="apeglm")

df <- as.data.frame(results_zm_t24Sig)

# Assuming 'df' is your data frame and 'column_name' is the name of the column
downr_df_sorted <- df[order(df$log2FoldChange), ]
head(downr_df_sorted)
downr_top <- head(downr_df_sorted, 20)

upr_df_sorted <- df[order(-df$log2FoldChange), ]
upr_top <- head(upr_df_sorted, 20)


df <- as.data.frame(downr_top)
df

return_results <- function(results, timepoint) {
  results <- subset(results, padj < 0.1)
  results_df <- as.data.frame(results)
  colnames_results <- lapply(colnames(results_df), function(x) paste0(x, timepoint))
  colnames(results_df) <- colnames_results
  return (results_df)
}

merge_dataframe <- function(first, second) {
  colnames(first)
  # first$gene_id <- rownames(first)
  second$gene_id <- rownames(second)
  head(second)
  merged_df <- merge(first, second, by.y = "gene_id", all.x = TRUE)
  return (merged_df)
}


merge_all_data <- function(main_df, other_dataframes,) {
  main_df$gene_id <- rownames(main_df)
  for (df in other_dataframes) {
    main_df <- merge_dataframe(main_df, df)
  }
  write.csv(main_df, file = glue("/Users/bsilke/bs-villungerlab/results/{filename}"))
}


# down_ndf <- ndf

library(tidyverse)

make_longdf_for_plot <- function(merged_df) {
  ndf <- data.frame(names<-merged_df$gene_id)
  ndf$n_16 <- merged_df$log2FoldChange_16
  ndf$n_20 <- merged_df$log2FoldChange_20
  ndf$n_24 <- merged_df$log2FoldChange
  ndf$n_36 <- merged_df$log2FoldChange_36
  ndf$n_48 <- merged_df$log2FoldChange_48
  ndf$symbol <- merged_df$symbol
  head(ndf)
  
  df_long <- ndf %>%
    pivot_longer(
      cols = starts_with("n_"), # Select columns that start with "n_"
      names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
      values_to = "log2foldchange" # The values will go into a new column "log2change"
    ) %>%
    mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
           Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric
  
  return (df_long)
}

# Plotting

plot_longdf <- function(long_df, treatment) {
  ggplot(down_df_long, aes(x = Timepoint, y = log2foldchange, group = symbol, color = symbol)) +
    geom_point() +
    geom_line() +
    labs(title = glue("{treatment} Treatment: downregulated genes"),
         x = "time(hrs)",
         y = expression(paste(log[2](x), 'fold change'))) +
    theme(plot.title = element_text(hjust = 0.5), # Center the title
          plot.title.position = "plot")
}

# geom_text(aes(label=symbol), hjust=0, vjust=0)
# theme(legend.position = "none") # Remove legend to avoid overcrowding if you have many genes

