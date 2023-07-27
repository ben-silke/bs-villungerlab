


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


return_results <- function(dds, coef, timepoint_extn, model='apeglm') {
  # Coef = "timepoint_t16_vs_t0", 
  # model = 'apeglm'
  results <- lfcShrink(dds, coef=coef, type=model)
  results <- add_annotations_to_results(results)
  return (results)
}
