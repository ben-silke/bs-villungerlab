


library(DESeq2)
library(DESeq2)
library("apeglm")
library(org.Hs.eg.db)
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)

setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')
dds_ZM <- ddseq_ZM

data_file = "~/bs-villungerlab/results/output_encode_1to6/ZM_unfiltered_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_ZM_t16 <- get_unfiltered_results(dds_ZM, "timepoint_t16_vs_t0", "_16")
  results_ZM_t20 <- get_unfiltered_results(dds_ZM, "timepoint_t20_vs_t0", "_20")
  results_ZM_t24 <- get_unfiltered_results(dds_ZM, "timepoint_t24_vs_t0", "_24")
  results_ZM_t36 <- get_unfiltered_results(dds_ZM, "timepoint_t36_vs_t0", "_36")
  results_ZM_t48 <- get_unfiltered_results(dds_ZM, "timepoint_t48_vs_t0", "_48")
  save(results_ZM_t48, results_ZM_t16, results_ZM_t20, results_ZM_t24, results_ZM_t36, file=data_file)
}

results_ZM_t16_df <- as.data.frame(results_ZM_t16)
results_ZM_t20_df <- as.data.frame(results_ZM_t20)
results_ZM_t24_df <- as.data.frame(results_ZM_t24)
results_ZM_t36_df <- as.data.frame(results_ZM_t36)
results_ZM_t48_df <- as.data.frame(results_ZM_t48)

all_df_merged_df <- merge_all_data(results_ZM_t48_df, results_ZM_t16_df, results_ZM_t20_df, results_ZM_t24_df, results_ZM_t36_df, 'results/output_encode/ZM/unfiltered_apeglm_ZM_data.csv', 'full_join')

file_increase <- "results/output_encode/ZM/iregulon_analysis/ZM_gene_signature_3n_increase.csv"
file_decrease <- "results/output_encode/ZM/iregulon_analysis/ZM_gene_signature_3n_increase.csv"

table_increase <- read.table(file_increase, sep=',', header=TRUE)
table_increase$Target.Gene
