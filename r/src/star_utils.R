library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")

source("r/src/pca_utils.R")
source("r/src/utils.R")

parse_star_filename <- function(filename, file_prefix=default_file_prefix) {
  pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
  matches <- str_match(filename, pattern)
  time <- str_extract(filename, "(?<=_)[^_]*(?=_r)")
  replicate <- str_extract(filename, "(?<=_r)\\d+")
  
  if (replicate <= 3) {
    batch = 1
  } else {
    batch = 2
  }
  
  name <- matches[1]
  
  return(c(name, time, replicate, batch))
}

create_star_dataframe <- function(treatment_name, data_directory, times, replicates_list) {
  replicates <- unlist(lapply(replicates_list, function(i) { paste0("r", i)}))
  # Create the files
  files <- list()
  counts <- list()
  for (time in times) {
    for (replicate in replicates) {
      file_path <- (file.path(data_directory, glue("{treatment_name}_{time}_{replicate}_fc")))
      files <- append(files, file_path)
    }
  }
  files <- unlist(files)
  
  print(files)
  print(file.exists(files))
  data_frame <- data.frame(files=files, stringsAsFactors = FALSE)
  parsed_values <- do.call("rbind", lapply(data_frame$files, parse_star_filename))
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  data_frame$batch <- paste0('b', parsed_values[, 4])
  return(data_frame)
}

create_star_dds <- function(treatment, data_directory, times, replicates_list, annotation_file, batch_correction=TRUE, trim_data=TRUE, run_DESeq=TRUE) {
  data_frame <- create_star_dataframe(treatment, data_directory, times, replicates_list)
  print(colnames(data_frame))
  # Check if all files exist
  print(data_frame$files)
  file.exists(data_frame$files)
  
  stopifnot(file.exists(data_frame$files))
  
  annotations <- read.delim(annotation_file)
  
  feature_counts_data <- featureCounts(data_frame$files, annot.ext=annotations)
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  
  
  # NEED TO WORK OUT WHAT IS COUNTS AND WHAT IS COLDATA?
  if (batch_correction) {
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ batch + timepoint)
    
  } else {
    dds <- DESeqDataSet(gene_se, design = ~ timepoint)
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ timepoint)
  }
  
  if (trim_data) {
    print(nrow(dds))
    keep <- rowSums(counts(dds)) >= 1
    dds <- dds[keep,]
    print(nrow(dds))
  }
  
  if (!run_DESeq) {
    return (dds)
  }
  
  # Run DESeq method on the data structure
  dds <- DESeq(dds)
  
  results <- results(dds)
  print(resultsNames(dds))
  return (dds)
}

create_htseq_dataframe <- function(treatment_name, data_directory, times, replicates_list) {
  replicates <- unlist(lapply(replicates_list, function(i) { paste0("r", i)}))
  # Create the files
  files <- list()
  counts <- list()
  for (time in times) {
    for (replicate in replicates) {
      file_path <- (file.path(data_directory, glue("htseq_{treatment_name}_{time}_{replicate}.counts")))
      files <- append(files, file_path)
    }
  }
  files <- unlist(files)
  
  print(files)
  print(file.exists(files))
  data_frame <- data.frame(files=files, stringsAsFactors = FALSE)
  vec <- lapply(data_frame$files, parse_star_filename)
  parsed_values <- do.call("rbind", vec)
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  data_frame$batch <- paste0('b', parsed_values[, 4])
  return(data_frame)
}


create_htseq_ddseq <- function(treatment_name, data_directory, times, replicates_list) {
  star_data <- create_htseq_dataframe(treatment_name, data_directory, times, replicates_list)
  df <- data.frame()
  for (file in star_data$files) {
    data <- read.table(file)
    print(file)
    pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
    matches <- str_match(file, pattern)
    names <- c('gene_id', matches[1])
    colnames(data) <- names
    df <- merge(df, data, all=TRUE)
  }

  # Remove the first column because it breaks and readd it
  merged_df <- df[,-2]
  file <- star_data$files[1]
  data <- read.table(file)
  pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
  matches <- str_match(file, pattern)
  names <- c('gene_id', matches[1])
  names[2]
  colnames(data) <- names
  merged_df <- merge(merged_df, data, by='gene_id')
  merged_df <- merged_df[-(1:5), ]

  # Fix the names
  rownames(merged_df) <- merged_df$gene_id
  merged_df <- merged_df[,-1]
  matrix <- as.matrix(merged_df)
  star_data$files = NULL

  stopifnot(dim(star_data)[1] == dim(matrix)[2])
  
  dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = star_data,
                              design= ~ batch + timepoint)
  
  dds <- DESeq(dds)
  return (dds)
}

load_all_htseq_data <- function(file_path) {
  file = file.path(file_path)
  data <- read.table(file)
  names <- list()
  for (name in colnames(data)) {
    parsed_name <- parse_star_filename(name)
    names <- append(names, parsed_name[1])
  }
  names
  
  data_frame <- data.frame(names=unlist(names), stringsAsFactors = FALSE)
  data_frame
  parsed_values <- do.call("rbind", lapply(data_frame$names, parse_star_filename))
  parsed_values
  
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  data_frame$batch <- paste0('b', parsed_values[, 4])
  colnames(data) <- data_frame$names
  
  start = dim(data)[1]-5
  end = dim(data)[1]-1+1
  df <- data[-(start:end), ]
  
  matrix <- as.matrix(data)
  dds <- DESeqDataSetFromMatrix(countData = matrix,
                                colData = data_frame,
                                design= ~ batch + timepoint)
  

  print(nrow(dds))
  keep <- rowSums(counts(dds)) >= 1
  dds <- dds[keep,]
  print(nrow(dds))

  dds <- DESeq(dds)
  return (dds)
}


