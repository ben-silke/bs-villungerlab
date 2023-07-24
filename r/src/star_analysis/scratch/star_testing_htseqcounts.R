library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")


source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")



#####
# ##### main
times = c(0, 16, 20, 24, 36, 48)
# times = c(0, 8, 12, 16, 24, 48)

treatment <- "ZM"
# data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_STAR'))
# data_directory = file.path('/Users/bsilke/bs-villungerlab/data/etop_star_output/output_STAR')
data_directory = file.path('/Volumes/bs_external/lab_villunger/htseq_counts/output_htseq_counts')
data_directory


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

# create_star_dataframe <- function(treatment_name, data_directory, times, replicates_list) {
star_data <- create_htseq_dataframe("ZM", data_directory, times, 1:6)
star_data



file <- star_data$files[1]
merged_data <- read.table(file)
df <- data.frame()
files
for (file in star_data$files) {
  data <- read.table(file)
  print(file)
  pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
  matches <- str_match(file, pattern)
  names <- c('gene_id', matches[1])
  colnames(data) <- names
  df <- merge(df, data, all=TRUE)
  rm(data)
}
test_df <- df
df <- test_df
merged_df <- df[,-2]

df
## add column zero 
file <- "/Volumes/bs_external/lab_villunger/htseq_counts/output_htseq_counts/htseq_ZM_0_r1.counts"
data <- read.table(file)
data
print(file)
pattern <- "([A-Z]*)_([0-9]*)_r([0-9]*)"
matches <- str_match(file, pattern)
names <- c('gene_id', matches[1])
names[2]
colnames(data) <- names

merged_df <- merge(merged_df, data, by='gene_id')

# rm(data)

# df.drop('B', axis=1, inplace=True)
old_merged_df <- merged_df
trimmed_df <- merged_df[-(1:5), ]

# merged_df <- df
trimmed_df

merged_df <- trimmed_df
rownames(merged_df) <- merged_df$gene_id
# merged_df <- merged_df[, -2]
merged_df <- merged_df[,-1]
matrix <- as.matrix(merged_df)
star_data$files = NULL
# star_data <- star_data[-1,]

dim(star_data)
dim(matrix)



dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = star_data,
                              design= ~ timepoint)


annotation <- get_annotation(dds)

# Where dds is a DESeqDataSet object
# pcaExplorer(dds=dds,
            # annotation=annotation)

dds <- DESeq(dds)


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c('batch', 'timepoint'))


mat <- assay(vsd)
mm <- model.matrix(~timepoint, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup=c('batch', 'timepoint'))

results <- lfcShrink(dds, coef="timepoint_t48_vs_t0", type="apeglm")
summary(results)

DESeq2::plotMA(results, ylim=c(-3,3))

results <- add_annotations_to_results(results)
selected_genes <- as.character(results$symbol)

EnhancedVolcano(results,
                lab = selected_genes,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8,8),
                ylab = expression(paste('-Log'[10],' adj P')),
                title = 'Differential expression for results_Nutl_t8_t0_shrunk_apeglm',
                pCutoff = 0.5,
                # at least double, or less than half
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)

