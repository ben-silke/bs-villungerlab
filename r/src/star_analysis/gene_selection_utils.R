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


source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_analysis/star_utils.R")


return_results <- function(dds, coef, timepoint_extn, model='apeglm') {
  # Coef = "timepoint_t16_vs_t0", 
  # model = 'apeglm'
  results <- lfcShrink(dds, coef=coef, type=model)
  results <- add_annotations_to_results(results)
  
  results <- subset(results, padj < 0.1)
  
  results_df <- as.data.frame(results)

  colnames_results <- lapply(colnames(results_df), function(x) paste0(x, timepoint_extn))
  colnames(results_df) <- colnames_results
  return (results_df)
}

merge_dataframe <- function(first, second) {
  second$gene_id <- rownames(second)
  print(head(second))
  print(head(first))
  merged_df <- merge(first, second, by.y = "gene_id", all.x = TRUE)
  return (merged_df)
}


merge_all_data <- function(main_df, one, two, three, four, filename) {
  main_df$gene_id <- rownames(main_df)
  
  merged_df <- merge_dataframe(main_df,one)
  merged_df <- merge_dataframe(merged_df,two)
  merged_df <- merge_dataframe(merged_df,three)
  merged_df <- merge_dataframe(merged_df,four)
  
  write.csv(merged_df, file = filename)
  return (merged_df)
}



make_longdf_for_plot <- function(merged_df, main_time) {
  
  # Ensure that all rows have a label
  if (main_time == 24) {
    merged_df$symbol <- merged_df$symbol_24
  } else if (main_time == 16) {
    merged_df$symbol <- merged_df$symbol_16
  } else {
    merged_df$symbol <- merged_df$symbol_48
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_24[is.na(merged_df$symbol)]
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_16[is.na(merged_df$symbol)]
  }
  
  merged_df$label <- merged_df$symbol
  print(merged_df$symbol)
  
  print(head(merged_df))
  merged_df$label[is.na(merged_df$label)] <- merged_df$gene_id[is.na(merged_df$label)]
  merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$gene_id[is.na(merged_df$symbol)]
  
  
  ndf <- data.frame(names<-merged_df$gene_id)
  ndf$n_16 <- merged_df$log2FoldChange_16
  ndf$n_20 <- merged_df$log2FoldChange_20
  ndf$n_24 <- merged_df$log2FoldChange_24
  ndf$n_36 <- merged_df$log2FoldChange_36
  ndf$n_48 <- merged_df$log2FoldChange_48
  ndf$symbol <- merged_df$label
  head(ndf)
  print(ndf)
  
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

plot_longdf <- function(long_df, plot_title) {
  p <- ggplot(long_df, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
    geom_point() +
    geom_line() +
    labs(title = plot_title,
         x = "time(hrs)",
         y = expression(paste(log[2](x), ' fold change'))) +
    theme(plot.title = element_text(hjust = 0.5), # Center the title
          plot.title.position = "plot")

  return (p)
}

# geom_text(aes(label=symbol), hjust=0, vjust=0)
# theme(legend.position = "none") # Remove legend to avoid overcrowding if you have many genes

