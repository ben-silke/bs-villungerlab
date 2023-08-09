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
dir.create('lab_work/qPCR/ZM_qPCR')


# /Users/bsilke/bs-villungerlab/lab_work/qPCR/ZM_qPCR/20230807-p53_time_series-ZM_results_validation.xlsx
EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/ZM_qPCR/"
FILE_NAME = "20230807-p53_time_series-ZM_results_validation.xlsx"

df_ZM = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')
df_ZM <- transform_data_values(df_ZM)
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")
full_df_ZM <- data.frame()

for (target in targets) {
  long_df_ZM <- create_long_df_for_sample(df_ZM, 'New', target)
  print(dim(long_df_ZM))
  full_df_ZM <- rbind(full_df_ZM, long_df_ZM)
}
long_df_ZM
colnames(full_df_ZM)
full_df_ZM

df_ZM_long_subset_targets_plot <- plot_ly(
  full_df_ZM,
  x = ~timepoint,
  y = ~avg_ddct2,
  mode = "markers+lines",
  hoverinfo = "text",
  # text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~target_name
) %>%
  layout(title="ZM Treatment - new : qPCR data")

df_ZM_long_subset_targets_plot
saveWidget(df_ZM_long_subset_targets_plot, "lab_work/qPCR/ZM_qPCR/ZM_all_new_pcr.html")



EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/Noc_qPCR/"
FILE_NAME = "20230808-p53_time_series-Noc_results_validation.xlsx"

df_Noc = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')
df_Noc <- transform_data_values(df_Noc)
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1")
View(df_Noc)
full_df_Noc <- data.frame()
for (target in targets) {
  long_df_Noc <- create_long_df_for_sample(df_Noc, 'plus', target)
  print(dim(long_df_Noc))
  full_df_Noc <- rbind(full_df_Noc, long_df_Noc)
}
long_df_Noc
colnames(full_df_Noc)
full_df_Noc

View(full_df_Noc)

df_Noc_long_subset_targets_plot <- plot_ly(
  full_df_Noc,
  x = ~timepoint,
  y = ~avg_ddct2,
  mode = "markers+lines",
  hoverinfo = "text",
  split= ~target_name
) %>%
  layout(title="Noc Treatment - new: qPCR data")

df_Noc_long_subset_targets_plot
saveWidget(df_Noc_long_subset_targets_plot, "lab_work/qPCR/Noc_qPCR/Noc_pseven_all_pcr.html")

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

new_df_noc <- save_html_for_treatment(df_Noc, 'new', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_new_all_pcr.html", "Nocodazole Treatment - new", "Noc")
p7_df_noc <- save_html_for_treatment(df_Noc, 'pseven', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_p7_all_pcr.html", "Nocodazole Treatment - p7", "Noc")
plus_df_noc <- save_html_for_treatment(df_Noc, 'plus', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_plus_all_pcr.html", "Nocodazole Treatment - +", "Noc")

names = colnames(new_df_noc)
names[2] = "avg_ddct2_new"
colnames(new_df_noc) = names
names[2] = "avg_ddct2_pseven"
colnames(p7_df_noc) = names
names[2] = "avg_ddct2_plus"
colnames(plus_df_noc) = names

complete_noc <- new_df_noc
# complete_noc$avg_ddct2_pseven <- p7_df_noc$avg_ddct2_pseven
complete_noc$avg_ddct2_plus <- plus_df_noc$avg_ddct2_plus
# complete_noc$avg <- (complete_noc$avg_ddct2_pseven + complete_noc$avg_ddct2_plus + complete_noc$avg_ddct2_new)/3
complete_noc$avg <- (complete_noc$avg_ddct2_plus + complete_noc$avg_ddct2_new)/2

# complete_noc <- subset(complete_noc)
complete_noc_avg <- plot_ly(
  complete_noc,
  x = ~timepoint,
  y = ~avg,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", timepoint, "<br>avg: ", avg, "<br>avg_ddct2_plus: ", avg_ddct2_plus, "<br>avg_ddct2_new: ", avg_ddct2_new),
  split= ~target_name
) %>%
  layout(title="Noc Average qPCR (new & +)")
# saveWidget(df_long_subset_targets_plot, save_location)
complete_noc_avg
saveWidget(complete_noc_avg, "lab_work/qPCR/Noc_qPCR/noc_average_newplus.html")

View(complete_noc)



######
new_df_zm <- save_html_for_treatment(df_ZM, 'New', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_new_all_pcr.html", "ZM Treatment - new")
p7_df_zm <- save_html_for_treatment(df_ZM, 'pseven', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_p7_all_pcr.html", "ZM Treatment - p7")
plus_df_zm <- save_html_for_treatment(df_ZM, 'plus', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_plus_all_pcr.html", "ZM Treatment - +")
names = colnames(new_df_zm)
names[2] = "avg_ddct2_new"
colnames(new_df_zm) = names
names[2] = "avg_ddct2_pseven"
colnames(p7_df_zm) = names
names[2] = "avg_ddct2_plus"
colnames(plus_df_zm) = names

complete_zm <- new_df_zm
# complete_zm$avg_ddct2_pseven <- p7_df_zm$avg_ddct2_pseven
complete_zm$avg_ddct2_plus <- plus_df_zm$avg_ddct2_plus
# complete_zm$avg <- (complete_zm$avg_ddct2_pseven + complete_zm$avg_ddct2_plus + complete_zm$avg_ddct2_new)/3
complete_zm$avg <- (complete_zm$avg_ddct2_plus + complete_zm$avg_ddct2_new)/2

# complete_zm <- subset(complete_zm)
complete_zm_avg <- plot_ly(
  complete_zm,
  x = ~timepoint,
  y = ~avg,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", timepoint, "<br>avg: ", avg, "<br>avg_ddct2_plus: ", avg_ddct2_plus, "<br>avg_ddct2_new: ", avg_ddct2_new),
  split= ~target_name
) %>%
  layout(title="zm Average qPCR (new and +)")
# saveWidget(df_long_subset_targets_plot, save_location)
complete_zm_avg
saveWidget(complete_zm_avg, "lab_work/qPCR/ZM_qPCR/ZM_average_newplus.html")

View(complete_zm)


