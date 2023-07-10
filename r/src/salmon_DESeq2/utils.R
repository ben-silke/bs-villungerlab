library(glue)
library(stringr)

library("pcaExplorer")
library(tximeta)
library(DESeq2)
library(org.Hs.eg.db)


default_file_prefix <- "salmon_quant_"

parse_filename <- function(filename, file_prefix=default_file_prefix) {
  
  name <- str_extract(filename, "(?<=salmon_quant_)[^_]*")
  time <- str_extract(filename, "(?<=_)[^_]*(?=_r)")
  replicate <- str_extract(filename, "(?<=_r)\\d+")
  
  if (replicate <= 3) {
    batch = 1
  } else {
    batch = 2
  }
  
  name <- glue("{name}_t{time}_r{replicate}")
  
  return(c(name, time, replicate, batch))
}

create_treatment_data <- function(treatment_name, data_directory, times, file_prefix, replicates_list) {
  replicates <- unlist(lapply(replicates_list, function(i) { paste0("r", i)}))
  # Create the files
  files <- list()
  for (time in times) {
    for (replicate in replicates) {
      file_path <- (file.path(data_directory, glue::glue("{file_prefix}_{treatment_name}_{time}_{replicate}"), 'quant.sf'))
      files <- append(files, file_path)
    }
  }
  files <- unlist(files)
  data_frame <- data.frame(files=files, stringsAsFactors = FALSE)
  parsed_values <- do.call("rbind", lapply(data_frame$files, parse_filename))
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  data_frame$batch <- paste0('b', parsed_values[, 4])
  return(data_frame)
}

create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE) {
  data_frame <- create_treatment_data(treatment, salmon_data_directory, times, file_prefix, replicates_list)
  print(colnames(data_frame))
  # Check if all files exist
  stopifnot(file.exists(data_frame$files))
  summarized_experiement <- tximeta(data_frame)
  gene_se <- summarizeToGene(summarized_experiement)
  
  if (batch_correction) {
    dds <- DESeqDataSet(gene_se, design = ~ batch + timepoint)
  } else {
    dds <- DESeqDataSet(gene_se, design = ~ timepoint)
  }
  
  if (trim_data) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]  
  }
  
  # Run DESeq method on the data structure
  dds <- DESeq(dds)
  
  results <- results(dds)
  print(resultsNames(dds))
  return (dds)
}