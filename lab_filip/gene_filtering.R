
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

# Post Filtering _ split by increase and decrease
########


subset_increase <- subset(merged_df,log2FoldChange.shrunk_24 > 0 )
subset_decrease <- subset(merged_df,log2FoldChange.shrunk_24 < 0 )

noquote_gene_list_increase <- noquote(subset_increase$gene_id)
noquote_gene_list_decrease <- noquote(subset_decrease$gene_id)

noquote_gene_list_increase <- c('gene_name', noquote_gene_list_increase)
noquote_gene_list_decrease <- c('gene_name', noquote_gene_list_decrease)

write(noquote_gene_list_increase, file = "results/lab_filip/increase_gene_list_padj005.txt")
write(noquote_gene_list_decrease, file = "results/lab_filip/decrease_gene_list_padj005.txt")

#### Generate long Df
######
colnames(merged_df)

ndf <- data.frame(names<-merged_df$gene_id)
ndf$n_12 <- merged_df$log2FoldChange.shrunk_12
ndf$n_24 <- merged_df$log2FoldChange.shrunk_24
ndf$n_30 <- merged_df$log2FoldChange.shrunk_30
ndf$symbol <- merged_df$gene_id
head(ndf)
print(ndf)

df_long_time <- ndf %>%
  pivot_longer(
    cols = starts_with("n_"), # Select columns that start with "n_"
    names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
    values_to = "log2foldchange" # The values will go into a new column "log2change"
  ) %>%
  mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
         Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric


## correct for padj
padj_df <- data.frame(names<-merged_df$gene_id)
padj_df$padj_12 <- merged_df$padj_12
padj_df$padj_24 <- merged_df$padj_24
padj_df$padj_30 <- merged_df$padj_30
padj_df$symbol <- merged_df$gene_id

df_long_padj <- padj_df %>%
  pivot_longer(
    cols = starts_with("padj_"), # Select columns that start with "n_"
    names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
    values_to = "padj" # The values will go into a new column "padj"
  ) %>%
  mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
         Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric

print(df_long_padj)
df_long <- merge(df_long_time, df_long_padj)
df_long$names....merged_df.gene_id <- NULL
colnames(df_long)

# headers <- colnames(df_long)
# headers[1] ='symbol'
# colnames(df_long_subset_targets) <- headers

df_long_subset_targets_plot <- plot_ly(
  df_long,
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


##########
