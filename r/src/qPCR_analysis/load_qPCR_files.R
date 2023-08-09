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

source("~/bs-villungerlab/r/src/qPCR_analysis/qPCR_utils.R")
# 
dir.create('lab_work/qPCR/NOC_qPCR')


# /Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/20230808-p53_time_series-Noc_results_validation.xlsx
EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/NOC_qPCR/"
FILE_NAME = "20230808-p53_time_series-Noc_results_validation.xlsx"

df = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')





df <- transform_data_values(df)
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")

full_df <- data.frame()
for (target in targets) {
  long_df <- create_long_df_for_sample(df, 'new', target)
  full_df <- merge(full_df, long_df, all = TRUE)
}

colnames(full_df)

df_long_subset_targets_plot <- plot_ly(
  full_df,
  x = ~timepoint,
  y = ~avg_ddct2,
  mode = "markers+lines",
  hoverinfo = "text",
  # text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~target_name
) %>%
  layout(title="Noc Treatment: qPCR data")

df_long_subset_targets_plot
saveWidget(df_long_subset_targets_plot, "lab_work/qPCR/NOC_qPCR/Noc_all_pcr.html")
