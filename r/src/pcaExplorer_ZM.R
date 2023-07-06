library("pcaExplorer")



# Importing salmon output for analysis with deSeq2
library(glue)
user <- 'bsilke'

dir <- glue("/Users/{user}/bs-villungerlab")
setwd(dir)
getwd()
salmon_data_directory = file.path(dir, glue('data/organised/{treatment}/output_salmon'))
times = c(0,8,12,16,20,24,36,48)

create_treatment_data <- function(treatment_name, data_directory, times, file_prefix) {
  files <- lapply(times, function(time) { file.path(data_directory, glue::glue("{file_prefix}_{treatment_name}_{time}"), 'quant.sf')})
  timepoints <- lapply(times, function(time) { paste0("t", time)})
  names <- lapply(times, function(time) { paste0(glue("{treatment_name}_"), time)})
  replicates <- lapply(1:6, function(i) { paste0("r", i)})
  metadata <- data.frame(files=files, names=names, timepoint=timepoints, replicate=replicates, stringsAsFactors = FALSE)
  return(list(metadata=metadata, files=files,timepoints=timepoints,names=names,replicates=replicates))
}

ZM_treatment_data <- create_treatment_data('ZM', salmon_data_directory,c(0,16,20,24,36,48), "salmon_quant") 
ZM_data_frame <- ZM_treatment_data$metadata
ZM_data_frame

#activate tximeta and build the SummarizeExperiment (se) object
library(tximeta)
ZM_summarized_experiment <- tximeta(ZM_data_frame)

gse <- summarizeToGene(ZM_summarized_experiment)
library(DESeq2)
dds <- DESeqDataSet(gse, design = ~ ZM_data_frame$timepoint)
dds


# Where dds is a DESeqDataSet object
pcaExplorer(dds=dds)