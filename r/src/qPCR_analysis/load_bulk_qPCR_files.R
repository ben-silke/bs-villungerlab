# Ok it is better to make these reusable, 
# but i think for now its faster for me to just create files where i need them
library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library("fgsea")


# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_utils.R")
# 
dir.create('lab_work/qPCR/fulloutput')
EXTENSION = "lab_work/raw_data/"
file_one = "qPCR_nocOld_zmOld_zm-19_noc5.xlsx"
# /Users/bsilke/bs-villungerlab/lab_work/raw_data/qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx
file_two = "qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx"
# /Users/bsilke/bs-villungerlab/lab_work/raw_data/qPCR_zmNew_nocNew_nocPlus_zmPlus.xlsx
file_three = "qPCR_zmNew_nocNew_nocPlus_zmPlus.xlsx"

# df_one = read_excel(file.path(EXTENSION, file_one), sheet='0')
# df_two = read_excel(file.path(EXTENSION, file_two), sheet='0')
# df_three = read_excel(file.path(EXTENSION, file_three), sheet='0')

# View(df_one)
# 
# colnames(df_one)
# 
# df <- df_one %>% dplyr::select('Target', 'Sample', 'adj_cq')
# cols <- colnames(df)
# cols[3] = 'Cq'
# colnames(df) = cols
# colnames(df)
# View(df)
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1")

save_html_for_treatment <- function(file, biological_replicate, targets, save_location, title, treatment, extension="lab_work/raw_data/") {
  df <- read_excel(file.path(EXTENSION, file), sheet='0')
  df <- df %>% dplyr::select('Target', 'Sample', 'adj_cq')
  cols <- colnames(df)
  cols[3] = 'Cq'
  colnames(df) = cols
  
  full_df <- data.frame()
  for (target in targets) {
    long_df <- create_longdf(df, biological_replicate, target)
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

create_longdf <- function(df, sample_name, target_name, control="GAPDH") {
  
  # SETUP control dataframe
  control_df <- subset(df, Target==control)
  sample_control_df <- subset(control_df, grepl(sample_name, Sample))
  sample_name
  control_df
  sample_control_df
  
  averaged_control_df <- sample_control_df %>%
    group_by(Sample) %>%
    summarise(
      avg_Cq_control = mean(Cq, na.rm = TRUE),
    )
  
  target_df <- subset(df, Target==target_name)
  df
  target_df
  sample_target_df <- subset(target_df, grepl(sample_name, Sample))
  target_df
  sample_target_df
  
  target_colnames = colnames(sample_target_df)
  target_colnames
  target_colnames[1] = 'sample_target'
  target_colnames[3] = 'target_Cq'
  colnames(sample_target_df) = target_colnames
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
# qPCR_nocOld_zmOld_zm-19_noc5.xlsx
df_one <- save_html_for_treatment("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_pcr.html", "Nocodazole Treatment - old", "Noc")
df_two <- save_html_for_treatment("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_pcr.html", "ZM Treatment - old", "ZM")
df_three <- save_html_for_treatment("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_all_pcr.html", "ZM Treatment - zm19", "ZM")
df_four <- save_html_for_treatment("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_pcr.html", "Nocodazole Treatment - zm5", "Noc")

# qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx
df_5 <- save_html_for_treatment("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_v2_pcr.html", "Nocodazole Treatment - old", "Noc")
df_6 <- save_html_for_treatment("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_v2_pcr.html", "ZM Treatment - old", "ZM")
df_7 <- save_html_for_treatment("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_v2_pcr.html", "ZM Treatment - zm19", "ZM")
df_8 <- save_html_for_treatment("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_v2_pcr.html", "Nocodazole Treatment - zm5", "Noc")


file_three = "qPCR_zmNew_nocNew_nocPlus_zmPlus.xlsx"
df_9 <- save_html_for_treatment(file_three, 'noc_plus', targets, "lab_work/qPCR/fulloutput/Noc_plus_pcr.html", "Nocodazole Treatment - plus", "Noc")
df_10 <- save_html_for_treatment(file_three, 'zm_plus', targets, "lab_work/qPCR/fulloutput/zm_plus_pcr.html", "ZM Treatment - plus", "ZM")
df_11<- save_html_for_treatment(file_three, 'zm_new', targets, "lab_work/qPCR/fulloutput/zm_new_pcr.html", "ZM Treatment - new", "ZM")
df_12 <- save_html_for_treatment(file_three, 'noc_new', targets, "lab_work/qPCR/fulloutput/noc_new_pcr.html", "Nocodazole Treatment - new", "Noc")
