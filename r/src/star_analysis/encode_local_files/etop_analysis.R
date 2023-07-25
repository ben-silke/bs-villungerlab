library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)

source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")
source("r/src/star_analysis/star_utils.R")

load('~/bs-villungerlab/results/output_encode_1to6/Etop_star_data.RData')

results_Etop_t12 <- lfcShrink(dds_Etop, coef="timepoint_t12_vs_t0", type="apeglm")
results_Etop_t8 <- lfcShrink(dds_Etop, coef="timepoint_t8_vs_t0", type="apeglm")
results_Etop_t16 <- lfcShrink(dds_Etop, coef="timepoint_t16_vs_t0", type="apeglm")
results_Etop_t24 <- lfcShrink(dds_Etop, coef="timepoint_t24_vs_t0", type="apeglm")
results_Etop_t48 <- lfcShrink(dds_Etop, coef="timepoint_t48_vs_t0", type="apeglm")

# return_results <- function(results, timepoint) {

# merge_dataframe <- function(first, second) {

# merge_all_data <- function(main_df, other_dataframes,) {


# make_longdf_for_plot <- function(merged_df) {

# plot_longdf <- function(long_df, treatment) {
