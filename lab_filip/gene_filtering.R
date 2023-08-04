
library(readxl)
library(plotly)
library(htmlwidgets)
# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")


# Alter where the data is.
EXTENSION = "lab_filip/transcript_deseq2_data"
dir.create("results/lab_filip")

files = c(
  "DESeq2_1.18.1.NoA_50nM_12h_VS_Control.alpha05.xlsx",
  "DESeq2_1.18.1.NoA_50nM_24h_VS_Control.alpha05.xlsx",
  "DESeq2_1.18.1.NoA_50nM_30h_VS_Control.alpha05.xlsx"
)



df_12 = read_excel(file.path(EXTENSION, files[1]))
df_24 = read_excel(file.path(EXTENSION, files[2]))
df_30 = read_excel(file.path(EXTENSION, files[3]))


append_timepoint_to_df_col <- function(df, timepoint_extn) {
  # If you have already run this operation, just return the dataframe
  for (colname in colnames(df)) {
    if (timepoint_extn %in% colname) {
      return (df)
    }
  }
  colnames_results <- lapply(colnames(df), function(x) paste0(x, timepoint_extn))
  colnames_results[1] <- 'gene_id'
  colnames(df) <- colnames_results
  return (df)
}

df_12 <- append_timepoint_to_df_col(df_12, '_12')
df_24 <- append_timepoint_to_df_col(df_24, '_24')
df_30 <- append_timepoint_to_df_col(df_30, '_30')

View(df_30)

merged_df <- merge(df_12, df_24, by.y = 'gene_id', all = TRUE)
merged_df <- merge(merged_df, df_30, by.y = 'gene_id', all = TRUE)

View(merged_df)

merged_df
dim(merged_df)
subset_merged_df <- subset(merged_df, padj_12 < 0.05)
subset_merged_df <- subset(subset_merged_df, padj_24 < 0.05)
subset_merged_df <- subset(subset_merged_df, padj_30 < 0.05)

View(subset_merged_df)

# Post Filtering
noquote_gene_list <- noquote(subset_merged_df$gene_id)
noquote_gene_list <- c('gene_name', noquote_gene_list)
write(noquote_gene_list, file = "results/lab_filip/gene_list_padj005.txt")



df_long_subset_targets <- generate_complete_long_df(subset_merged_df, 24)
View(df_long_subset_targets)

headers <- colnames(df_long_subset_targets)
headers[1] ='symbol'
colnames(df_long_subset_targets) <- headers


df_long_subset_targets_plot <- plot_ly(
  df_long_subset_targets,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
) 
# %>%
  # layout(title="Transcription Dataset")

df_long_subset_targets_plot
saveWidget(df_long_subset_targets_plot, "results/lab_filip/transcript_data.html")

