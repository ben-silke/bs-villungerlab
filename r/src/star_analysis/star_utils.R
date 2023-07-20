library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")


parse_star_filename <- function(filename, file_prefix=default_file_prefix) {
  pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
  matches <- str_match(filename, pattern)
  # name <- str_extract(filename, "(?<=.counts)[^_]*")
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