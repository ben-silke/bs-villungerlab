

# Ok it is better to make these reusable, 
# but i think for now its faster for me to just create files where i need them
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

source("~/bs-villungerlab/r/src/qPCR_analysis/gene_selection_utils.R")


# /Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/20230808-p53_time_series-Noc_results_validation.xlsx
EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/"
FILE_NAME = "20230808-p53_time_series-Noc_results_validation.xlsx"

df = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')

# View(df)
df$Content = NULL

# df$Sample <- gsub("+","p", df$Sample)

df$Sample <- gsub(" 0h","_0h", df$Sample)
df$Sample <- gsub(" 16h","_16h", df$Sample)
df$Sample <- gsub(" 20h","_20h", df$Sample)
df$Sample <- gsub(" 24h","_24h", df$Sample)
df$Sample <- gsub(" 36h","_36h", df$Sample)
df$Sample <- gsub(" 48h","_48h", df$Sample)

df$Cq <- as.numeric(as.character(df$Cq))

control_df <- subset(df, Target=="GAPDH")
# View(control_df)
colnames(control_df)
sample_new_df <- subset(control_df, grepl('new', Sample))
plus_new_df <- subset(control_df, grepl('+', Sample))
p7_new_df <- subset(control_df, grepl('p7', Sample))
# View(control_df)

df
groupd_df <- df %>%
  group_by(Target)
groupd_df

avg_control_df <- control_df %>%
  group_by(Sample) %>%
  summarise(avg_Cq = mean(Cq, na.rm = TRUE))

averaged_df <- df %>%
  group_by(Sample) %>%
  summarise(
    avg_Cq = mean(Cq, na.rm = TRUE),
  )


create_graph_for_sample <- function(sample_name, target_name, control="GAPDH") {
  control="GAPDH"
  sample_name="new"
  
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
  # sample_control_df
  # sample_control_df$Target = control
  # sample_control_df
  
  avg_colnames = colnames(sample_control_df)
  avg_colnames
  avg_colnames[3] = 'control_target'
  avg_colnames[5] = 'control_Cq'
  colnames(sample_control_df) = avg_colnames
  sample_control_df$Fluor = NULL
  sample_control_df$Well = NULL
  sample_control_df
  # Now for the specific treatment
  target_name = 'BMF'
  # 
  
  target_df <- subset(df, Target==target_name)
  sample_target_df <- subset(target_df, grepl(sample_name, Sample))
  target_df
  sample_target_df
  
  # averaged_target_df <- sample_target_df %>%
  #   group_by(Sample) %>%
  #   summarise(
  #     avg_Cq_target = mean(Cq, na.rm = TRUE),
  #   )
  # averaged_target_df
  # averaged_target_df$Target = target_name
  # averaged_target_df
  
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
  
  avg_ddCt2_df <- next_merged_df %>%
    group_by(Sample) %>%
    summarise(
      avg_ddct2 = mean(ddct2, na.rm = TRUE),
    )
  
  final_merged_df <- left_join(next_merged_df, avg_dCt_df)
  final_merged_df
  # merged_df$dCt <- merged_df$avg_Cq_target - merged_df$avg_Cq_control
  # merged_df
  
  long_df <- final_merged_df %>%
    group_by(Sample) %>%
    summarise(
      avg_ddct2 = mean(ddct2, na.rm = TRUE),
    )
  long_df$target_name <-   target_name
  
  long_df
  return (long_df)
}

df_plot <- plot_ly(
  long_df,
  x = ~Sample,
  y = ~avg_ddct2,
  mode = "markers+lines",
  hoverinfo = "text",
  # text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  # split= ~symbol
) %>%
  layout(title="ZM Treatment: Primers")
df_plot
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")
