
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")

# Load the dplyr package
library(dplyr)

source("r/src/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

file <- file.path("data/allgc_ZM_fc")


df <- read.csv(file, sep="\t", row.names="Geneid", skip=1)
df
df$Chr <- NULL
df$Start <- NULL
df$End <- NULL
df$Strand <- NULL
df$Length <- NULL
df
matrix <- as.matrix(df)

names <- colnames(matrix)
short_names <- names[6:length(names)]


# table <- data_frame(sample<- )
# coldata <- table()


# Extract the time and replicate information
time <- (str_extract(names, "(?<=ZM_)[0-9]+"))
replicate <- (str_extract(names, "(?<=_r)[0-9]+"))


# Combine the time and replicate information into a data frame (table)
table <- data.frame(file=names, time = as.numeric(time), replicate = as.numeric(replicate))
table$batch <- ifelse(table$replicate <= 3, "1", "2")

table <- table %>%
  mutate_all(~as.character(.))

dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = table,
                              design = ~ batch + time)

dds
print(nrow(dds))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
print(nrow(dds))
dds <- estimateSizeFactors(dds)

# counts(dds) <- counts(dds) + 1
counts(dds)[counts(dds) == 0] <- 1L

# counts(dds)[counts(dds) == 0] <- 0.005

dds <- DESeq(dds)

results <- results(dds)
print(resultsNames(dds))

plotMA(results)
#############
fix_colnames <- function(input_string) {
  zm_substring <- sub(".*gc_(ZM)_.*", "\\1", input_string)
  num_substring <- sub(".*gc_ZM_(\\d+)_.*", "\\1", input_string)
  r_substring <- sub(".*gc_ZM_\\d+_(r\\d+).*", "\\1", input_string)
  
  return (glue("{zm_substring}_{num_substring}_{r_substring}"))
}

new_names <- lapply(names, fix_colnames)
names
new_names

nm <- unlist(new_names)
nm

colnames(matrix) <- nm

matrix
