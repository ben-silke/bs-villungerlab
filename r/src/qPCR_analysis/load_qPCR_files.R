

# Ok it is better to make these reusable, 
# but i think for now its faster for me to just create files where i need them
library(readxl)
library(plotly)
library(htmlwidgets)
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

control_df <- subset(df, Target=="GAPDH")



df_plot <- plot_ly(
  df,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
) %>%
  layout(title="ZM Treatment: Primers")


zm_df <- fix_labels(all_df_merged_df)


"FOXM1" %in% zm_df$symbol 
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")
subset_targets <- subset(zm_df, symbol %in% targets)
subset_targets

df_long_subset_targets <- generate_complete_long_df(subset_targets, 24)


df_long_subset_targets_plot
saveWidget(df_long_subset_targets_plot, "results/output_encode/ZM/ZM_all_primers.html")

