# TODO: make this non-specific and see if you can make it like click

library(glue)
library(stringr)

library("pcaExplorer")
library(tximeta)
library(DESeq2)
library(org.Hs.eg.db)

user <- 'bsilke'
dir <- glue("/Users/{user}/bs-villungerlab")
setwd(dir)
getwd()

treatment <- "ZM"
salmon_data_directory = file.path(dir, glue('data/organised/{treatment}/output_salmon'))
times = c(0,8,12,16,20,24,36,48)
times <- c(0,16,24,36,48)

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
  # timepoints <- unlist(c(lapply(times, function(time) { paste0("t", time)})))
  # names <- unlist(lapply(times, function(time) { paste0(glue("{treatment_name}_"), time)}))
  data_frame <- data.frame(files=files, stringsAsFactors = FALSE)
  parsed_values <- do.call("rbind", lapply(data_frame$files, parse_filename))
  data_frame$names <- parsed_values[, 1]
  data_frame$timepoint <- paste0("t", parsed_values[, 2])
  data_frame$replicate <- paste0("r", parsed_values[, 3])
  data_frame$batch <- paste0('b', parsed_values[, 4])
  return(data_frame)
}

create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE) {
  data_frame <- create_treatment_data(treatment, salmon_data_directory, times, files_prefix, replicates_list)
  print(colnames(data_frame))
  # Check if all files exist
  stopifnot(file.exists(data_frame$files))
  summarized_experiement <- tximeta(data_frame)
  gene_se <- summarizeToGene(summarized_experiement)
  
  if (batch_correction) {
    dds <- DESeqDataSet(gse_one, design = ~ batch + timepoint)
  } else {
    dds <- DESeqDataSet(gse_one, design = ~ timepoint)
  }
  
  if (trim_data) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]  
  }
  
  dds <- DESeq(dds)
  results <- results(dds)
  print(resultsNames(dds))
  return (dds)
}

ZM_data_frame_exp_one <- create_treatment_data('ZM', salmon_data_directory, times, "salmon_quant", 1:6)
colnames(ZM_data_frame_exp_one)
# ZM_data_frame_exp_two <- create_treatment_data('ZM', salmon_data_directory, times, "salmon_quant", 4:6) 

#Check that all files exist
stopifnot(file.exists(ZM_data_frame_exp_one$files))

#activate tximeta and build the SummarizeExperiment (se) object
ZM_summarized_experiment_one <- tximeta(ZM_data_frame_exp_one)
# ZM_summarized_experiment_two <- tximeta(ZM_data_frame_exp_two)

gse_one <- summarizeToGene(ZM_summarized_experiment_one)
# gse_two <- summarizeToGene(ZM_summarized_experiment_two)

dds_one <- DESeqDataSet(gse_one, design = ~ batch + timepoint)
# dds_two <- DESeqDataSet(gse_two, design = ~ 1)

# Get the annotation
# library(org.Hs.eg.db)
# genenames <- mapIds(org.Hs.eg.db,keys = rownames(dds_one),column = "SYMBOL",keytype="ENSEMBL")
# # genenames
# annotation <- data.frame(gene_name=genenames,
#                          row.names=rownames(dds_one),
#                          stringsAsFactors = FALSE)
# 
# head(annotation) 

keep <- rowSums(counts(dds_one)) >= 10
dds_one <- dds_one[keep,]

# dds_one$condition <- factor(dds$condition, levels= ZM_data_frame_exp_one$replicate)
dds <- dds_one
# Analysis with DESeq2
dds <- DESeq(dds)

results <- results(dds)
results

resultsNames(dds)



# zero_sixteen <- results(dds, contrast=c("timepoint", "t0", "t16"))

ordered_results <- results[order(results$pvalue),]
ordered_results
summary(ordered_results)

sum(ordered_results$padj < 0.05, na.rm=TRUE)

results_alpha05 = results(dds, alpha=0.05)
summary(results_alpha05)


plotMA(results, ylim=c(-2,2))


# plotCounts(dds, gene=which.min(results$padj), intgroup="timepoint")


##### results_t16_v_t0 
results_t16_v_t0 <- results(dds, name="timepoint_t16_vs_t0")
results_t16_v_t0

results_t16_v_t0_LFC <- lfcShrink(dds, coef="timepoint_t16_vs_t0", type="apeglm")
results_t16_v_t0_LFC

results_t16_v_t0 <- results_t16_v_t0[order(results_t16_v_t0$pvalue),]
results_t16_v_t0
summary(results_t16_v_t0)
sum(results_t16_v_t0$padj < 0.05, na.rm=TRUE)

plotMA(results_t16_v_t0, ylim=c(-2,2))
plotMA(results_t16_v_t0_LFC, ylim=c(-2,2))

## Using Ashr
results_t16_v_t0_LFC <- lfcShrink(dds, coef="timepoint_t16_vs_t0", type="apeglm")
results_t16_v_t0_LFC

results_t16_v_t0 <- results_t16_v_t0[order(results_t16_v_t0$pvalue),]
results_t16_v_t0
summary(results_t16_v_t0)
sum(results_t16_v_t0$padj < 0.05, na.rm=TRUE)

# plotMA(results_t16_v_t0, ylim=c(-2,2))
# plotMA(results_t16_v_t0_LFC, ylim=c(-2,2))


# plotCounts(dds, gene=which.min(results$padj), intgroup="timepoint")


## Count data transformations
## blind should be used if we are looking at individual batches, 
# but if we are looking at everything then we want to consider the difference which is explainable by the batch effect which we have identified in the PCA
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


## heatmap for count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


## likelihood ratio test:
# wald test or lRT test
dds <- DESeq(dds, test="LRT", reduced=~batch)
# this accounts for condition as well
res <- results(dds)
res
