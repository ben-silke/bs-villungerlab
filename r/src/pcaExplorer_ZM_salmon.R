
library(glue)
library(stringr)

library("pcaExplorer")
library(tximeta)
library(DESeq2)

user <- 'bsilke'

dir <- glue("/Users/{user}/bs-villungerlab")
setwd(dir)
getwd()

treatment <- "ZM"
salmon_data_directory = file.path(dir, glue('data/organised/{treatment}/output_salmon'))
times = c(0,8,12,16,20,24,36,48)
times <- c(0,16)

default_file_prefix <- "salmon_quant_"

parse_filename <- function(filename, file_prefix=default_file_prefix) {

  name <- str_extract(filename, "(?<=salmon_quant_)[^_]*")
  time <- str_extract(filename, "(?<=_)[^_]*(?=_r)")
  replicate <- str_extract(filename, "(?<=_r)\\d+")
  
  name <- glue("{name}_t{time}_r{replicate}")
  
  return(c(name, time, replicate))
}

create_treatment_data <- function(treatment_name, data_directory, times, file_prefix) {
  replicates <- unlist(lapply(1:6, function(i) { paste0("r", i)}))
  # Create the files
  files <- list()
  for (time in times) {
    for (replicate in replicates) {
      file_path <- (file.path(data_directory, glue::glue("{file_prefix}_{treatment_name}_{time}_{replicate}"), 'quant.sf'))
      files <- append(files, file_path)
    }
  }
  files <- unlist(files)
  # timepoints <- unlist(c(lapply(times, function(time) { paste0("t", time)})))
  # names <- unlist(lapply(times, function(time) { paste0(glue("{treatment_name}_"), time)}))
  data_frame <- data.frame(files=files, stringsAsFactors = FALSE)
  parsed_values <- do.call("rbind", lapply(data_frame$files, parse_filename))
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  
  return(data_frame)
}

ZM_data_frame <- create_treatment_data('ZM', salmon_data_directory, times, "salmon_quant") 
ZM_data_frame

#Check that all files exist
stopifnot(file.exists(ZM_data_frame$files))

#activate tximeta and build the SummarizeExperiment (se) object
ZM_summarized_experiment <- tximeta(ZM_data_frame)

gse <- summarizeToGene(ZM_summarized_experiment)

dds <- DESeqDataSet(gse, design = ~ 1)
dds

# Where dds is a DESeqDataSet object
pcaExplorer(dds=dds)
