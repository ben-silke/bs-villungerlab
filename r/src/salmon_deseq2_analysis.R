# Importing salmon output for analysis with deSeq2
library(glue)
user <- 'bsilke'

dir <- glue("/Users/{user}/bs-villungerlab")
setwd(dir)
getwd()



treatment <- 'ZM'
data_directory = file.path(dir, glue('data/organised/{treatment}/output_salmon'))
times = c(0,8,12,16,20,24,36,48)
files <- lapply(times, function(time) { file.path(data_directory, glue::glue("salmon_quant_{treatment}_{time}"), 'quant.sf')})

files
names <- lapply(times, function(i) { paste0("ZM_", i)})
#creating the timepoint variable
timepoints <- lapply(times, function(i) { paste0("t", i)})
# timepoint <- c("t0","t8","t12","t16", "t24", "t48")
#creating replicate variable
replicates <- lapply(1:6, function(i) { paste0("r", i)})
replicates

#building the coldata dataframe
coldata <- data.frame(files = files, names= names, timepoint=timepoints, replicate = replicates, stringsAsFactors=FALSE)
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
