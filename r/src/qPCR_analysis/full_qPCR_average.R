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
source("~/bs-villungerlab/r/src/qPCR_analysis/load_bulk_qPCR_files.R")


new_df_noc <- save_html_for_treatment(df_Noc, 'new', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_new_all_pcr.html", "Nocodazole Treatment - new", "Noc")
p7_df_noc <- save_html_for_treatment(df_Noc, 'pseven', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_p7_all_pcr.html", "Nocodazole Treatment - p7", "Noc")
plus_df_noc <- save_html_for_treatment(df_Noc, 'plus', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/Noc_qPCR/Noc_plus_all_pcr.html", "Nocodazole Treatment - +", "Noc")

new_df_zm <- save_html_for_treatment(df_ZM, 'New', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_new_all_pcr.html", "ZM Treatment - new")
p7_df_zm <- save_html_for_treatment(df_ZM, 'pseven', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_p7_all_pcr.html", "ZM Treatment - p7")
plus_df_zm <- save_html_for_treatment(df_ZM, 'plus', c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ANKRD1"), "lab_work/qPCR/ZM_qPCR/ZM_plus_all_pcr.html", "ZM Treatment - +")

# qPCR_nocOld_zmOld_zm-19_noc5.xlsx
noc_old <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_pcr.html", "Nocodazole Treatment - old", "Noc")
zm_old <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_pcr.html", "ZM Treatment - old", "ZM")
zm_nineteen <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_all_pcr.html", "ZM Treatment - zm19", "ZM")
noc_pfive <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_pcr.html", "Nocodazole Treatment - zm5", "Noc")

# qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx
noc_old_v2 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_v2_pcr.html", "Nocodazole Treatment - old", "Noc")
zm_old_v2 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_v2_pcr.html", "ZM Treatment - old", "ZM")
zm_nineteen_v2 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_v2_pcr.html", "ZM Treatment - zm19", "ZM")
noc_pfive_v2 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_v2_pcr.html", "Nocodazole Treatment - zm5", "Noc")

file_three = "qPCR_zmNew_nocNew_nocPlus_zmPlus.xlsx"
noc_plus <- save_html_file(file_three, 'noc_plus', targets, "lab_work/qPCR/fulloutput/Noc_plus_pcr.html", "Nocodazole Treatment - plus", "Noc")
zm_plus <- save_html_file(file_three, 'zm_plus', targets, "lab_work/qPCR/fulloutput/zm_plus_pcr.html", "ZM Treatment - plus", "ZM")
zm_new <- save_html_file(file_three, 'zm_new', targets, "lab_work/qPCR/fulloutput/zm_new_pcr.html", "ZM Treatment - new", "ZM")
noc_new <- save_html_file(file_three, 'noc_new', targets, "lab_work/qPCR/fulloutput/noc_new_pcr.html", "Nocodazole Treatment - new", "Noc")


# Merge All files into one dataframe, take averages,
# plot as one timeseries.