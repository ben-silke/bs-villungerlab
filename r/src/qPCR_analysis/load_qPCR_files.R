# Ok it is better to make these reusable, 
# but i think for now its faster for me to just create files where i need them
library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_utils.R")
# 
dir.create('lab_work/qPCR/ZM_qPCR')
EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/ZM_qPCR/"
FILE_NAME = "20230807-p53_time_series-ZM_results_validation.xlsx"
df_ZM = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')
df_ZM <- transform_data_values(df_ZM)

EXTENSION = "/Users/bsilke/bs-villungerlab/lab_work/qPCR/Noc_qPCR/"
FILE_NAME = "20230808-p53_time_series-Noc_results_validation.xlsx"
df_Noc = read_excel(file.path(EXTENSION, FILE_NAME), sheet='raw_data')
df_Noc <- transform_data_values(df_Noc)


new_df_noc <- save_html_for_treatment(df_Noc, 'new', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_new_all_pcr.html", "Nocodazole Treatment - new", "Noc")
new_df_noc
p7_df_noc <- save_html_for_treatment(df_Noc, 'pseven', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_p7_all_pcr.html", "Nocodazole Treatment - p7", "Noc")
p7_df_noc
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
complete_noc <- complete_noc %>%
  rowwise() %>%
  mutate(sd_ddct2 = sd(c(avg_ddct2_plus, avg_ddct2_new)))

complete_noc <- complete_noc %>%
  rowwise() %>%
  mutate(mean_ddct2 = mean(c(avg_ddct2_plus, avg_ddct2_new)))

View(complete_noc)

complete_noc_avg <- plot_ly(
  complete_noc,
  x = ~timepoint,
  y = ~mean_ddct2,
  # error_y = list(array = ~sd_ddct2, color = "red", thickness = 0.5, width = 2),
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", timepoint, "<br>mean: ", mean_ddct2, "<br>sd: ", sd_ddct2, "<br>avg_ddct2_plus: ", avg_ddct2_plus, "<br>avg_ddct2_new: ", avg_ddct2_new),
  split= ~target_name
) %>%
  layout(title="Noc Average qPCR (all)")
  # saveWidget(df_long_subset_targets_plot, save_location)
complete_noc_avg
saveWidget(complete_noc_avg, "lab_work/qPCR/Noc_qPCR/noc_all_average.html")

View(complete_noc)

complet

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
