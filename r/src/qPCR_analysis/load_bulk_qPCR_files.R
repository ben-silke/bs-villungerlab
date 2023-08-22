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
source("~/bs-villungerlab/r/src/qPCR_analysis/qPCR_utils.R")
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


# qPCR_nocOld_zmOld_zm-19_noc5.xlsx
df_one <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_pcr.html", "Nocodazole Treatment - old", "Noc")
df_two <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_pcr.html", "ZM Treatment - old", "ZM")
df_three <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_all_pcr.html", "ZM Treatment - zm19", "ZM")
df_four <- save_html_file("qPCR_nocOld_zmOld_zm-19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_pcr.html", "Nocodazole Treatment - zm5", "Noc")

# qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx
df_5 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_old', targets, "lab_work/qPCR/fulloutput/Noc_old_v2_pcr.html", "Nocodazole Treatment - old", "Noc")
df_6 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_old', targets, "lab_work/qPCR/fulloutput/zm_old_v2_pcr.html", "ZM Treatment - old", "ZM")
df_7 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'zm_nineteen', targets, "lab_work/qPCR/fulloutput/zm_19_v2_pcr.html", "ZM Treatment - zm19", "ZM")
df_8 <- save_html_file("qPCR_v2_nocOld_zmOld_zm19_noc5.xlsx", 'noc_pfive', targets, "lab_work/qPCR/fulloutput/noc_5_v2_pcr.html", "Nocodazole Treatment - zm5", "Noc")


file_three = "qPCR_zmNew_nocNew_nocPlus_zmPlus.xlsx"
df_9 <- save_html_file(file_three, 'noc_plus', targets, "lab_work/qPCR/fulloutput/Noc_plus_pcr.html", "Nocodazole Treatment - plus", "Noc")
df_10 <- save_html_file(file_three, 'zm_plus', targets, "lab_work/qPCR/fulloutput/zm_plus_pcr.html", "ZM Treatment - plus", "ZM")
df_11<- save_html_file(file_three, 'zm_new', targets, "lab_work/qPCR/fulloutput/zm_new_pcr.html", "ZM Treatment - new", "ZM")
df_12 <- save_html_file(file_three, 'noc_new', targets, "lab_work/qPCR/fulloutput/noc_new_pcr.html", "Nocodazole Treatment - new", "Noc")
