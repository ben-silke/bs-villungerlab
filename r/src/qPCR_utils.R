library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(dplyr)

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")

create_long_df_for_sample <- function(df, sample_name, target_name, control="GAPDH") {

  # SETUP control dataframe
  control_df <- subset(df, Target==control)
  sample_control_df <- subset(control_df, grepl(sample_name, Sample))
  control_df
  sample_control_df
  
  averaged_control_df <- sample_control_df %>%
    group_by(Sample) %>%
    summarise(
      avg_Cq_control = mean(Cq, na.rm = TRUE),
    )
  
  target_df <- subset(df, Target==target_name)
  sample_target_df <- subset(target_df, grepl(sample_name, Sample))
  target_df
  sample_target_df
  
  target_colnames = colnames(sample_target_df)
  target_colnames
  target_colnames[3] = 'sample_target'
  target_colnames[5] = 'target_Cq'
  colnames(sample_target_df) = target_colnames
  sample_target_df$Fluor = NULL
  sample_target_df$Well = NULL
  sample_target_df$SQ = NULL
  sample_target_df
  averaged_control_df
  
  # TODO: Handle target_Cq values which are zero, or less than a threshold better (where there is a clear error)
  # merged_df <- left_join(sample_target_df, sample_control_df, by = c('Sample'))
  merged_df <- merge(sample_target_df, averaged_control_df)
  merged_df
  merged_df$dCt <- merged_df$target_Cq - merged_df$avg_Cq_control
  merged_df
  
  avg_dCt_df <- merged_df %>%
    group_by(Sample) %>%
    summarise(
      avg_dCt_target = mean(dCt, na.rm = TRUE),
    )
  avg_dCt_df
  next_merged_df <- left_join(merged_df, avg_dCt_df)
  next_merged_df
  
  dCt_avg_t0 <- avg_dCt_df$avg_dCt_target[1]
  
  next_merged_df$ddCt <- next_merged_df$dCt - dCt_avg_t0
  next_merged_df$ddct2 <- 2^-(next_merged_df$ddCt)
  print(next_merged_df)
  
  avg_ddCt2_df <- next_merged_df %>%
    group_by(Sample) %>%
    summarise(
      avg_ddct2 = mean(ddct2, na.rm = TRUE),
    )
  
  final_merged_df <- left_join(next_merged_df, avg_dCt_df)
  print(final_merged_df)
  
  long_df <- final_merged_df %>%
    group_by(Sample) %>%
    summarise(
      avg_ddct2 = mean(ddct2, na.rm = TRUE),
    )
  long_df$target_name <- target_name
  
  long_df
  long_df$timepoint <- as.numeric(str_extract(long_df$Sample, "\\d+"))
  long_df
  return (long_df)
}


transform_data_values <- function(df) {
  print(df)
  colnames(df)
  df$Content = NULL
  df$SQ = NULL
  df$Cq = NULL
  df
  # df$Cq = df$adj_cq
  df$Sample <- gsub(" 0h","_0h", df$Sample)
  df$Sample <- gsub(" 16h","_16h", df$Sample)
  df$Sample <- gsub(" 20h","_20h", df$Sample)
  df$Sample <- gsub(" 24h","_24h", df$Sample)
  df$Sample <- gsub(" 36h","_36h", df$Sample)
  df$Sample <- gsub(" 48h","_48h", df$Sample)
  
  names <- colnames(df)
  names[5] <- "Cq"
  colnames(df) <- names
  df$Cq <- as.numeric(as.character(df$Cq))
  df
  return (df)
}


save_html_for_treatment <- function(df, biological_replicate, targets, save_location, title, treatment) {
  full_df <- data.frame()
  for (target in targets) {
    long_df <- create_long_df_for_sample(df, biological_replicate, target)
    full_df <- rbind(full_df, long_df)
  }
  df_long_subset_targets_plot <- plot_ly(
    full_df,
    x = ~timepoint,
    y = ~avg_ddct2,
    mode = "markers+lines",
    text = ~paste("T: ", timepoint, "<br>avg_ddct2: ", avg_ddct2),
    hoverinfo = "text",
    split= ~target_name
  ) %>%
    layout(title=title)
  saveWidget(df_long_subset_targets_plot, save_location)
  full_df$treatment = treatment
  return(full_df)
}
