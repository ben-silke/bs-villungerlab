# Importing salmon output for analysis with deSeq2
library(glue)
user <- 'bsilke'

dir <- glue("/Users/{user}/bs-villungerlab")
setwd(dir)
getwd()

treatment <- 'Nutl'
data_directory = file.path(dir, glue('data/organised/{treatment}/output_salmon'))
Nutl_t0 <- file.path(data_directory,"salmon_quant_Nutl_0","quant.sf")
Nutl_t8 <- file.path(data_directory,"salmon_quant_Nutl_8","quant.sf")
Nutl_t12 <- file.path(data_directory,"salmon_quant_Nutl_12","quant.sf")
Nutl_t16 <- file.path(data_directory,"salmon_quant_Nutl_16","quant.sf")
Nutl_t24 <- file.path(data_directory,"salmon_quant_Nutl_24","quant.sf")
Nutl_t48 <- file.path(data_directory,"salmon_quant_Nutl_48","quant.sf")

files <- c(Nutl_t0, Nutl_t8, Nutl_t12, Nutl_t16, Nutl24, Nutl48)
names <- c("Nutl_t0", "Nutl_t8", "Nutl_t12", "Nutl_t16", "Nutl24", "Nutl48")
#creating the timepoint variable
timepoint <- c("t0","t8","t12","t16", "t24", "t48")
#creating replicate variable
replicate <- c("r1")

#building the coldata dataframe
coldata <- data.frame(files = files, names= names, timepoint=timepoint, replicate = replicate, stringsAsFactors=FALSE)
coldata

#activate tximeta and build the SummarizeExperiment (se) object
library(tximeta)
se <- tximeta(coldata)

1#make the se (summarize experiment) ready for DGE analysis on Deseq2
gse <- summarizeToGene(se)
library(DESeq2)

dds <- DESeqDataSet(gse, design = ~ timepoint)


nrow(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)


# DESeq requires replicates
dds <- DESeq(dds)
results <- results(dds)
results
