library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")

source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

#####
# ##### main
times = c(0, 16, 20, 24, 36, 48)
times = c(0, 8, 12, 16, 24, 48)
treatment <- "Etop"
# data_directory = file.path('../../../Volumes/bs_external/lab_villunger', glue('organised/{treatment}/output_STAR'))
# data_directory = file.path('/Users/bsilke/bs-villungerlab/data/etop_star_output/output_STAR')
data_directory = file.path('/nobackup/lab_villunger/bsilke/organised/Etop/output_STAR')
data_directory

# Can we do this on the cluster ??
# annotation_file <- "/Users/bsilke/bs-villungerlab/data/hg38_RefSeq_22Apr2022.saf"

annotation_file <- "/nobackup/lab_villunger/bsilke/references/Homo_sapiens.GRCh38.109.gtf"

##### bam_files
# /Volumes/nobackup/lab_villunger/bsilke/organised/ZM/output_STAR/ZM_12_r5Aligned.out.sam

# files_feature_counts <- create_star_dds('ZM', data_directory, times, 1:6)
files_feature_counts <- create_star_dataframe(treatment, data_directory, times, 1)

# getwd()
T2T <- read.delim(annotation_file)
T2T
type(files_feature_counts$files)
# files_feature_counts$files
is.vector(files_feature_counts$files)

fc <- featureCounts(files_feature_counts$files, annot.ext=T2T)

dds <- DESeqDataSetFromMatrix(countData = fc$counts,
                              colData = files_feature_counts,
                              design = ~ condition)
dds


# Check dds works
dds <- DESeq(dds)
results <- results(dds)
print(resultsNames(dds))







