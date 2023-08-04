
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
dir.create("results/lab_filip")

EXTENSION = "lab_lukas/Omics/RNAseq/"
FILE_NAME = "top50_differential.xlsx"


df = read_excel(file.path(EXTENSION, FILE_NAME))




noquote_gene_list <- noquote(df$gene_id)
noquote_gene_list <- c('gene_name', noquote_gene_list)
write(noquote_gene_list, file = "results/lab_lukas/top_50_differential.txt")

