library(glue)
library(stringr)

library("pcaExplorer")
library(tximeta)
library(DESeq2)
library(org.Hs.eg.db)
library("RUVSeq")

source("r/src/pca_utils.R")


DEFAULT_TIMES_SHORT <- c(0, 8, 16, 20, 24, 48)
DEFAULT_TIMES_LONG <- c(0, 16, 20, 24, 36, 48)

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

create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE, run_DESeq=TRUE) {
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
    keep <- rowSums(counts(dds)) >= 1
    dds <- dds[keep,]  
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

add_annotations_to_results <- function(res) {
  ens.str <- substr(rownames(res), 1, 15)
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  return (res)
}

create_htseq_count_df <-
  function(treatment_name,
           data_directory,
           times,
           file_prefix,
           replicates_list) {
    replicates <-
      unlist(lapply(replicates_list, function(i) {
        paste0("r", i)
      }))
    # Create the files
    files <- list()
    for (time in times) {
      for (replicate in replicates) {
        file_path <-
          (file.path(
            data_directory,
            glue::glue("{file_prefix}_{treatment_name}_{time}_{replicate}"),
            glue(
              'htseq_count_{treatment_name}_{time}_{replicate}.counts'
            )
          ))
        files <- append(files, file_path)
      }
    }
    files <- unlist(files)
    data_frame <-
      data.frame(files = files, stringsAsFactors = FALSE)
    parsed_values <-
      do.call("rbind", lapply(data_frame$files, parse_filename))
    data_frame$names <- parsed_values[, 1]
    data_frame$timepoint <- paste0("t", parsed_values[, 2])
    data_frame$replicate <- paste0("r", parsed_values[, 3])
    data_frame$batch <- paste0('b', parsed_values[, 4])
    return(data_frame)
  }

create_dds_from_htseq <-
  function(treatment,
           ht_seq_data_directory,
           times,
           file_prefix,
           replicates_list,
           batch_correction = TRUE,
           trim_data = TRUE,
           run_DESeq = TRUE) {
    data_frame <-
      create_treatment_data(treatment,
                            ht_seq_data_directory,
                            times,
                            file_prefix,
                            replicates_list)
    print(colnames(data_frame))
    # Check if all files exist
    stopifnot(file.exists(data_frame$files))
    
    
    if (batch_correction) {
      ddsHTseq <- DESeqDataSetFromHTSeqCount(
        sampleTable = data_frame,
        directory = ht_seq_data_directory,
        design = ~ batch + timepoint
      )
    } else {
      ddsHTseq <- DESeqDataSetFromHTSeqCount(
        sampleTable = data_frame,
        directory = ht_seq_data_directory,
        design = ~ timepoint
      )
    }
    
    if (trim_data) {
      keep <- rowSums(counts(ddsHTseq)) >= 10
      ddsHTseq <- ddsHTseq[keep,]
    }
    
    if (!run_DESeq) {
      return (ddsHTseq)
    }
    
    # Run DESeq method on the data structure
    ddsHTseq <- DESeq(ddsHTseq)
    
    results <- results(ddsHTseq)
    print(resultsNames(ddsHTseq))
    return (ddsHTseq)
  }


account_for_batch_effect <- function(dds) {
  set <- newSeqExpressionSet(counts(dds))
  idx  <- rowSums(counts(set) > 5) >= 2
  set  <- set[idx, ]
  set <- betweenLaneNormalization(set, which="upper")
  not.sig <- rownames(res)[which(res$pvalue > .1)]
  empirical <- rownames(set)[ rownames(set) %in% not.sig ]
  set <- RUVg(set, empirical, k=2)
  print(pData(set))
  return (set)
}

