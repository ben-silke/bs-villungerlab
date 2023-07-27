library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")

source("r/src/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

#####
# ##### main
times = c(0, 16, 20, 24, 36, 48)
# times = c(0, 8, 12, 16, 24, 48)

treatment <- "ZM"
# data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_STAR'))
# data_directory = file.path('/Users/bsilke/bs-villungerlab/data/etop_star_output/output_STAR')
data_directory = file.path('data/zm_fc_counts/output_feature_counts')
data_directory

# create_star_dataframe <- function(treatment_name, data_directory, times, replicates_list) {
star_data <- create_star_dataframe("ZM", data_directory, times, 1:6)
star_data$files

matrix_list <- list()
for (file in star_data$files) {
  df <- read.csv(file, sep="\t", row.names="Geneid", skip=1)
  matrix <- as.matrix(df)
  matrix_list <- append(matrix_list, matrix)
}

all_data <- do.call(cbind, matrix_list)

head(all_data)

# time <- 0
# replicate <- "r1"
# file_path <- (file.path(data_directory, glue("{treatment}_{time}_{replicate}_fc")))
# 
# count_data <- read.csv(file_path, sep="\t", row.names="Geneid", skip=1)
# count_data
# 
# cts <- as.matrix(count_data)
# print(cts)
# files <- append(files, file_path)
# counts <- append(counts, cts)


######## 
# create_star_dds <- function(treatment, data_directory, times, replicates_list, annotation_file, batch_correction=TRUE, trim_data=TRUE, run_DESeq=TRUE) {
#   data_frame <- create_star_dataframe(treatment, data_directory, times, replicates_list)
#   print(colnames(data_frame))
#   # Check if all files exist
#   print(data_frame$files)
#   file.exists(data_frame$files)
#   
#   stopifnot(file.exists(data_frame$files))
#   
#   annotations <- read.delim(annotation_file)
#   
#   feature_counts_data <- featureCounts(data_frame$files, annot.ext=annotations)
#   
#   dds <- DESeqDataSetFromMatrix(countData = cts,
#                                 colData = coldata,
#                                 design = ~ condition)
#   
#   
#   
#   # NEED TO WORK OUT WHAT IS COUNTS AND WHAT IS COLDATA?
#   if (batch_correction) {
#     dds <- DESeqDataSetFromMatrix(countData = cts,
#                                   colData = coldata,
#                                   design = ~ batch + timepoint)
#     
#   } else {
#     dds <- DESeqDataSet(gene_se, design = ~ timepoint)
#     dds <- DESeqDataSetFromMatrix(countData = cts,
#                                   colData = coldata,
#                                   design = ~ timepoint)
#   }
#   
#   if (trim_data) {
#     print(nrow(dds))
#     keep <- rowSums(counts(dds)) >= 1
#     dds <- dds[keep,]
#     print(nrow(dds))
#   }
#   
#   if (!run_DESeq) {
#     return (dds)
#   }
#   
#   # Run DESeq method on the data structure
#   dds <- DESeq(dds)
#   
#   results <- results(dds)
#   print(resultsNames(dds))
#   return (dds)
# }




read_featureCounts <- function( path=".", pattern, reshape=TRUE, stats=FALSE){
  if(stats){
    if(missing(pattern))  pattern <- "\\.summary$"
    # second column  is always unique, so skip header and assign column names
    fc <- read_sample_files(path, pattern, col_names=c("status", "count"), skip=1)
    if(reshape){
      fc <- dplyr::filter(fc, count!=0) %>% tidyr::spread(status, count)
      #  add option to sort HCI samples as X1, X2, ..., X10, X11
      fc  <-  fc[ order_samples(fc$sample), ]
    }
  }
  else{
    if(missing(pattern))  pattern <- "\\.counts$"
    fc <- read_sample_files(path, pattern, col_names=c("geneid",	"chr", "start",	"end", "strand", "length", "count"), skip=2)
    if(reshape){
      fc  <- dplyr::select(fc, sample, geneid, count) %>% tidyr::spread(sample, count)
      fc  <-  fc[, c(1, order_samples(colnames(fc)[-1])+1) ]
    }
  }
  fc
}
