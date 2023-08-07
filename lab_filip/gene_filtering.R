
library(readxl)
library(plotly)
library(htmlwidgets)
library("docstring")

# Import necessary util functions and set the directory correctly.
setwd("~/bs-villungerlab/"). # If you have cloned the repo to your home folder, and the directory structure is the same this will work. otherwise, change
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")


# Alter where the data is if your directory structure is different.
EXTENSION = "lab_filip/transcript_deseq2_data"
dir.create("results/lab_filip/12_hours")
dir.create("results/lab_filip/12_hours")

files = c(
  "DESeq2_1.18.1.NoA_50nM_12h_VS_Control.alpha05.xlsx",
  "DESeq2_1.18.1.NoA_50nM_24h_VS_Control.alpha05.xlsx",
  "DESeq2_1.18.1.NoA_50nM_30h_VS_Control.alpha05.xlsx"
)

# Read files into data frame format. If you are dealing with ddseq objects these should be results objects, which you can then use.
df_12 = read_excel(file.path(EXTENSION, files[1]))
df_24 = read_excel(file.path(EXTENSION, files[2]))
df_30 = read_excel(file.path(EXTENSION, files[3]))


append_timepoint_to_df_col <- function(df, timepoint_extn) {
  #' Helper function to append the string, timepoint_extn to each header/ colname in the dataframe
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

merged_df <- merge(df_12, df_24, by.y = 'gene_id', all = TRUE)
merged_df <- merge(merged_df, df_30, by.y = 'gene_id', all = TRUE)

merged_df
dim(merged_df)
subset_merged_df <- subset(merged_df, padj_12 < 0.05)
subset_merged_df <- subset(subset_merged_df, padj_24 < 0.05)
subset_merged_df <- subset(subset_merged_df, padj_30 < 0.05)

dim(subset_merged_df)
# greater abs(lfc) than 1 for which time point? 12hrs
subset_merged_df <- subset(subset_merged_df, abs(log2FoldChange.shrunk_12) > 1)

dim(subset_merged_df)

# Post Filtering _ split by increase and decrease
########
subset_increase <- subset(subset_merged_df,log2FoldChange.shrunk_24 > 0 )
subset_decrease <- subset(subset_merged_df,log2FoldChange.shrunk_24 < 0 )

dim(subset_increase)
dim(subset_decrease)


noquote_gene_list_increase <- noquote(subset_increase$gene_id)
noquote_gene_list_decrease <- noquote(subset_decrease$gene_id)

noquote_gene_list_increase <- c('gene_name', noquote_gene_list_increase)
noquote_gene_list_decrease <- c('gene_name', noquote_gene_list_decrease)

write(noquote_gene_list_increase, file = "results/lab_filip/12_hours/increase_gene_list_padj005.txt")
write(noquote_gene_list_decrease, file = "results/lab_filip/12_hours/decrease_gene_list_padj005.txt")

#### Generate long Df
######
colnames(subset_merged_df)

create_long_df <- function(merged_df) {
  #' Helper function to create a long dataframe from the flat dataframe with padj and lfc values to allow timeseries.
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
  return (df_long)
}


df_long_increase <- create_long_df(subset_increase)

df_long_subset_targets_plot_increase <- plot_ly(
  df_long_increase,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
) %>%
  layout(title="12 hours differentiallly expressed genes:  |LFC|>1 && adjp<0.05")

df_long_subset_targets_plot_increase
saveWidget(df_long_subset_targets_plot_increase, "results/lab_filip/12_hours/df_long_subset_targets_plot_increase.html")

##########


df_long_decrease <- create_long_df(subset_decrease)

df_long_subset_targets_plot_decrease <- plot_ly(
  df_long_decrease,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
) %>%
  layout(title="12 hours differentiallly expressed genes:  |LFC|>1 && adjp<0.05")

df_long_subset_targets_plot_decrease
saveWidget(df_long_subset_targets_plot_decrease, "results/lab_filip/12_hours/df_long_subset_targets_plot_decrease.html")

