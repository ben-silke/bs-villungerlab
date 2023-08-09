library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(dplyr)

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_analysis/qPCR_utils.R")


create_long_df_for_sample <- function(sample_name, target_name, control="GAPDH") {
  
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
  
  avg_colnames = colnames(sample_control_df)
  avg_colnames
  avg_colnames[3] = 'control_target'
  avg_colnames[5] = 'control_Cq'
  colnames(sample_control_df) = avg_colnames
  
  sample_control_df$Fluor = NULL
  sample_control_df$Well = NULL
  sample_control_df
  target_df <- subset(df, Target==target_name)
  sample_target_df <- subset(target_df, grepl(sample_name, Sample))
  target_df
  sample_target_df
  
  target_colnames = colnames(sample_target_df)
  target_colnames[3] = 'sample_target'
  target_colnames[5] = 'target_Cq'
  colnames(sample_target_df) = target_colnames
  sample_target_df$Fluor = NULL
  sample_target_df$Well = NULL
  
  sample_target_df
  averaged_control_df
  
  # merged_df <- left_join(sample_target_df, sample_control_df, by = c('Sample'))
  merged_df <- merge(sample_target_df, averaged_control_df)
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
  
  next_merged_df$ddCt <- next_merged_df$dCt - next_merged_df$avg_dCt_target
  next_merged_df$ddct2 <- 2^next_merged_df$ddCt
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
  long_df$target_name <-target_name
  
  long_df
  long_df$timepoint <- as.numeric(str_extract(long_df$Sample, "\\d+"))
  long_df
  return (long_df)
}