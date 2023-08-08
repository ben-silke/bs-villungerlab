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
library(docstring)


source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")


return_results <- function(dds, coef, timepoint_extn, model='apeglm') {
  #' Function to shrink results with DESeq2 from a dds object according to a specific model.
  #' #' returns results as a dataframe with the timepoint appended to the headers/ colnames.
  #' Coef = "timepoint_t16_vs_t0", 
  #' model = 'apeglm'
  results <- lfcShrink(dds, coef=coef, type=model)
  results <- add_annotations_to_results(results)
  
  results <- subset(results, padj < 0.1)
  
  results_df <- as.data.frame(results)

  colnames_results <- lapply(colnames(results_df), function(x) paste0(x, timepoint_extn))
  colnames(results_df) <- colnames_results
  return (results_df)
}

merge_dataframe <- function(first, second, join_type='none', by_y="gene_id") {
  #' Function to merge two dataframes with specific join condition.
  second$gene_id <- rownames(second)
  print(dim(second))
  print(dim(first))
  print(typeof(first))
  print(typeof(second))
  if (join_type == 'full_join') {
    merged_df <- merge(first, second, by.y = by_y, all = TRUE)
  } else {
    merged_df <- merge(first, second, by.y = by_y, all.x = TRUE)
  }
  return (merged_df)
}


# This is probably the wrong way to do this,
# TODO: implement better handling of multiple dataframes, or with an list/ array format
merge_all_data <- function(main_df, one, two=FALSE, three=FALSE, four=FALSE, filename, join_type='') {
  main_df$gene_id <- rownames(main_df)
  print(rownames(main_df))

  merged_df <- merge_dataframe(main_df, one, join_type)

  if (two & is.data.frame(two)) { merged_df <- merge_dataframe(merged_df, two, join_type) }  
  if (three & is.data.frame(three)) { merged_df <- merge_dataframe(merged_df, three, join_type) }  
  if (four & is.data.frame(four)) { merged_df <- merge_dataframe(merged_df, four, join_type) }  
  
  write.csv(merged_df, file = filename)
  return (merged_df)
}



make_longdf_for_plot <- function(merged_df, main_time) {
  # Ensure that all rows have a label/ symbol
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
          plot.title.position = "plot",
          legend.position = "bottom")
  return (p)
}

get_unfiltered_results <- function(dds, coef, timepoint_extn, model='apeglm') {
  results <- lfcShrink(dds, coef=coef, type=model)
  results <- add_annotations_to_results(results)
  results_df <- as.data.frame(results)
  colnames_results <- lapply(colnames(results_df), function(x) paste0(x, timepoint_extn))
  colnames(results_df) <- colnames_results
  return (results_df)
}

append_timepoint <- function(results_df, timepoint_extn) {
  colnames_results <- lapply(colnames(results_df), function(x) paste0(x, timepoint_extn))
  colnames(results_df) <- colnames_results
  return (results_df)
}


fix_labels <- function(df) {
  df$symbol <- df$symbol_48
  df$symbol[is.na(df$symbol)] <- df$symbol_16[is.na(df$symbol)]
  df$symbol[is.na(df$symbol)] <- df$symbol_24[is.na(df$symbol)]
  df$symbol[is.na(df$symbol)] <- df$gene_id[is.na(df$symbol)]
  return (df)
}


generate_complete_long_df <- function(merged_df, main_time) {
  #' generates long dataframe from flat dataframe with padj and lfc change
  # main_time = 24
  # Ensure that all rows have a label
  if (main_time == 24) {
    merged_df$symbol <- merged_df$symbol_24
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_48[is.na(merged_df$symbol)]
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_16[is.na(merged_df$symbol)]
  } else if (main_time == 16) {
    merged_df$symbol <- merged_df$symbol_16
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_48[is.na(merged_df$symbol)]
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_24[is.na(merged_df$symbol)]
  } else if (main_time == 'special') {
    merged_df$symbol = merged_df$gene_id
  } else {
    merged_df$symbol <- merged_df$symbol_48
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_24[is.na(merged_df$symbol)]
    merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_16[is.na(merged_df$symbol)]
  }
  merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$gene_id[is.na(merged_df$symbol)]
  
  
  merged_df$label <- merged_df$symbol
  print(merged_df$symbol)
  
  print(head(merged_df))
  merged_df$label[is.na(merged_df$label)] <- merged_df$gene_id[is.na(merged_df$label)]
  merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$gene_id[is.na(merged_df$symbol)]
  
  ## correct for log2foldchange
  ndf <- data.frame(names<-merged_df$gene_id)
  ndf$n_16 <- merged_df$log2FoldChange_16
  ndf$n_20 <- merged_df$log2FoldChange_20
  ndf$n_24 <- merged_df$log2FoldChange_24
  ndf$n_36 <- merged_df$log2FoldChange_36
  ndf$n_48 <- merged_df$log2FoldChange_48
  ndf$symbol <- merged_df$label
  head(ndf)
  print(ndf)
  
  df_long_time <- ndf %>%
    pivot_longer(
      cols = starts_with("n_"), # Select columns that start with "n_"
      names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
      values_to = "log2foldchange" # The values will go into a new column "log2change"
    ) %>%
    mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
           Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric
  
  
  ## correct for padj
  padj_df <- data.frame(names<-merged_df$gene_id)
  padj_df$padj_16 <- merged_df$padj_16
  padj_df$padj_20 <- merged_df$padj_20
  padj_df$padj_24 <- merged_df$padj_24
  padj_df$padj_36 <- merged_df$padj_36
  padj_df$padj_48 <- merged_df$padj_48
  padj_df$symbol <- merged_df$label
  head(padj_df)
  print(padj_df)
  
  df_long_padj <- padj_df %>%
    pivot_longer(
      cols = starts_with("padj_"), # Select columns that start with "n_"
      names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
      values_to = "padj" # The values will go into a new column "padj"
    ) %>%
    mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
           Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric
  
  print(df_long_padj)
  df_long <- merge(df_long_time, df_long_padj)
  df_long  
  return (df_long)
}