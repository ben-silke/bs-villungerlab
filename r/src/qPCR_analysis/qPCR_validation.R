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

dir.create("lab_work/qPCR/validation")
EXTENSION = "lab_work/qPCR/raw_data"
# noc_data_frame <- read_file("noc_validation_output_data.csv", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_v2_pcr.html", "Nocodazole Treatment - old", "Noc", extension=EXTENSION)

# /Users/bsilke/bs-villungerlab/lab_work/qPCR/raw_data/noc_validation_outputdata.xlsx
# /Users/bsilke/bs-villungerlab/lab_work/qPCR/raw_data/zm_target_validation.xlsx
noc_df <- read_excel(file.path(EXTENSION, "noc_target_validation.xlsx"), sheet='0')
zm_df <- read_excel(file.path(EXTENSION, "zm_target_validation.xlsx"), sheet='0')

targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "AURKB", "CDKN1A")


ZM_RPE1 <- read_file("zm_target_validation.xlsx", 'RPE1', targets, "lab_work/qPCR/validation/RPE1_qPCR_validation.html", "ZM Treatment - RPE1", "ZM")
ZM_A549 <- read_file("zm_target_validation.xlsx", 'A549', targets, "lab_work/qPCR/validation/A549_qPCR_validation.html", "ZM Treatment - A549", "ZM")
ZM_nalm6_wt <- read_file("zm_target_validation.xlsx", 'Nalm6_WT', targets, "lab_work/qPCR/validation/nalm6_wt_qPCR_validation.html", "ZM Treatment - Nalm6 WT", "ZM")
ZM_nalm6_hc9 <- read_file("zm_target_validation.xlsx", 'Nalm6_HC9', targets, "lab_work/qPCR/validation/nalm6_hc9_qPCR_validation.html", "ZM Treatment - Nalm6 HC9", "ZM")


full_df <- rbind(ZM_RPE1, ZM_A549, ZM_nalm6_wt, ZM_nalm6_hc9)

full_df
write.csv(full_df,'lab_work/qPCR/validation/full_zm_validation.csv')


targets = c("ARID5B", 'BMF', "CDC25A", "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "AURKB", "CDKN1A")


Noc_RPE1 <- read_file("Noc_target_validation.xlsx", 'RPE1', targets, "lab_work/qPCR/validation/RPE1_qPCR_validation.html", "Noc Treatment - RPE1", "Noc")
Noc_A549 <- read_file("Noc_target_validation.xlsx", 'A549', targets, "lab_work/qPCR/validation/A549_qPCR_validation.html", "Noc Treatment - A549", "Noc")
Noc_nalm6_wt <- read_file("Noc_target_validation.xlsx", 'Nalm6_WT', targets, "lab_work/qPCR/validation/nalm6_wt_qPCR_validation.html", "Noc Treatment - Nalm6 WT", "Noc")
Noc_nalm6_hc9 <- read_file("Noc_target_validation.xlsx", 'Nalm6_HC9', targets, "lab_work/qPCR/validation/nalm6_hc9_qPCR_validation.html", "Noc Treatment - Nalm6 HC9", "Noc")


full_df <- rbind(Noc_RPE1, Noc_A549, Noc_nalm6_wt, Noc_nalm6_hc9)

full_df
write.csv(full_df,'lab_work/qPCR/validation/full_Noc_validation.csv')

long_df
# save_bar_chart

View(ZM_RPE1) 


process_data <- function(df, sample_name, target_name, control="GAPDH") {
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
  averaged_control_df
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
  
  # avg_ddCt2_df <- next_merged_df %>%
  #   group_by(Sample) %>%
  #   summarise(
  #     avg_ddct2 = mean(ddct2, na.rm = TRUE),
  #     var2_ddct2 = var(ddct2,na.rm = TRUE)^2
  #   )
  
  final_merged_df <- left_join(next_merged_df, avg_dCt_df)
  print(final_merged_df)
  final_merged_df$target_name <- target_name
  print(final_merged_df)
  
  return(final_merged_df)
  
  # long_df <- final_merged_df %>%
  #   group_by(Sample) %>%
  #   summarise(
  #     avg_ddct2 = mean(ddct2, na.rm = TRUE),
  #     sd_ddct2 = sd(ddct2, na.rm = TRUE),
  #     var_ddct2 = var(ddct2, na.rm = TRUE),
  #     var2_ddct2 = var(ddct2,na.rm = TRUE)^2
  #   )
  # long_df$target_name <- target_name
  # 
  # long_df
  # return (long_df)
}





# df <- zm_df
# cell_type="A549"
# target_name ="ZMAT3"
# control="GAPDH"
read_file <- function(file, cell_type, targets, save_location, title, treatment, extension="lab_work/raw_data/") {
  df <- read_excel(file.path(EXTENSION, file), sheet='0')
  df <- df %>% dplyr::select('Target', 'Sample', 'adj_cq')
  cols <- colnames(df)
  cols[3] = 'Cq'
  colnames(df) = cols
  
  full_df <- data.frame()
  target= target_name
  for (target in targets) {
    long_df <- process_data(df, cell_type, target)
    long_df$treatment = substr(long_df$Sample, start=1, stop=nchar(long_df$Sample)-2)
    long_df
    full_df <- rbind(full_df, long_df)
  }
  long_df
  full_df
  
  return(full_df)
}

df_new
compute_shared_sd <- function(mu1, sigma1, mu2, sigma2) {
  n1=4
  n2=4
  combined_var <- ((n1 - 1)*sigma1^2 + (n2 - 1)*sigma2^2 + 
                     n1*n2*(mu1 - mu2)^2 / (n1 + n2)) / (n1 + n2 - 1)
  combined_sd <- sqrt(combined_var)
  return (combinded_sd)
}




