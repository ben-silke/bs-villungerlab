library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)


load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')

dds_zm <- ddseq_ZM

results_zm_t24 <- lfcShrink(dds_zm, coef="timepoint_t24_vs_t0", type="apeglm")
results_zm_t24Sig <- subset(results_zm_t24, padj < 0.1)
results_zm_t24Sig <- add_annotations_to_results(results_zm_t24Sig)

results_zm_t16 <- lfcShrink(dds_zm, coef="timepoint_t16_vs_t0", type="apeglm")
results_zm_t20 <- lfcShrink(dds_zm, coef="timepoint_t20_vs_t0", type="apeglm")
results_zm_t36 <- lfcShrink(dds_zm, coef="timepoint_t36_vs_t0", type="apeglm")
results_zm_t48 <- lfcShrink(dds_zm, coef="timepoint_t48_vs_t0", type="apeglm")

df <- as.data.frame(results_zm_t24Sig)

# Assuming 'df' is your data frame and 'column_name' is the name of the column
downr_df_sorted <- df[order(df$log2FoldChange), ]
head(downr_df_sorted)
downr_top <- head(downr_df_sorted, 20)

upr_df_sorted <- df[order(-df$log2FoldChange), ]
upr_top <- head(upr_df_sorted, 20)


df <- as.data.frame(downr_top)
df


results_zm_t16Sig <- subset(results_zm_t16, padj < 0.1)
res_16_df <- as.data.frame(results_zm_t16Sig)
colnames(res_16_df)
typeof(colnames(res_16_df))
# colnames_res16 <- lapply(colnames(res_16_df), paste0(str, "_16"))
colnames_res16 <- lapply(colnames(res_16_df), function(x) paste0(x, "_16"))
colnames_res16
colnames(res_16_df) <- colnames_res16
res_16_df

results_zm_t20Sig <- subset(results_zm_t20, padj < 0.1)
res_20_df <- as.data.frame(results_zm_t20Sig)
colnames(res_20_df)
typeof(colnames(res_20_df))
# colnames_res20 <- lapply(colnames(res_20_df), paste0(str, "_20"))
colnames_res20 <- lapply(colnames(res_20_df), function(x) paste0(x, "_20"))
colnames_res20
colnames(res_20_df) <- colnames_res20
res_20_df


# results_zm_t24 <- lfcShrink(dds_zm, coef="timepoint_t24_vs_t0", type="apeglm")
# results_zm_t24Sig <- subset(results_zm_t24, padj < 0.1)

results_zm_t36Sig <- subset(results_zm_t36, padj < 0.1)
res_36_df <- as.data.frame(results_zm_t36Sig)
colnames(res_36_df)
typeof(colnames(res_36_df))
# colnames_res36 <- lapply(colnames(res_36_df), paste0(str, "_36"))
colnames_res36 <- lapply(colnames(res_36_df), function(x) paste0(x, "_36"))
colnames_res36
colnames(res_36_df) <- colnames_res36
res_36_df

results_zm_t48Sig <- subset(results_zm_t48, padj < 0.1)
res_48_df <- as.data.frame(results_zm_t48Sig)
colnames(res_48_df)
typeof(colnames(res_48_df))
# colnames_res48 <- lapply(colnames(res_48_df), paste0(str, "_48"))
colnames_res48 <- lapply(colnames(res_48_df), function(x) paste0(x, "_48"))
colnames_res48
colnames(res_48_df) <- colnames_res48
res_48_df

colnames(df)
df$gene_id <- rownames(df)
res_16_df$gene_id <- rownames(res_16_df)
head(df)
merged_df <- df_merged <- merge(df, res_16_df, by.y = "gene_id", all.x = TRUE)
merged_df
head(merged_df)

merge_dataframe <- function(first, second) {
  colnames(first)
  # first$gene_id <- rownames(first)
  second$gene_id <- rownames(second)
  head(second)
  merged_df <- merge(first, second, by.y = "gene_id", all.x = TRUE)
  return (merged_df)
}

merged_df <- merge_dataframe(merged_df,res_20_df)
merged_df <- merge_dataframe(merged_df,res_36_df)
merged_df <- merge_dataframe(merged_df,res_48_df)

colnames(merged_df)

# write.csv(merged_df, file = "/Users/bsilke/bs-villungerlab/results/merged_up_regulated_zm.csv")
write.csv(merged_df, file = "/Users/bsilke/bs-villungerlab/results/down_regulated_zm.csv")

ndf <- data.frame(names<-merged_df$gene_id)
ndf$n_16 <- merged_df$log2FoldChange_16
ndf$n_20 <- merged_df$log2FoldChange_20
ndf$n_24 <- merged_df$log2FoldChange
ndf$n_36 <- merged_df$log2FoldChange_36
ndf$n_48 <- merged_df$log2FoldChange_48
ndf$symbol <- merged_df$symbol
head(ndf)


# down_ndf <- ndf

library(tidyverse)

df_long <- ndf %>%
  pivot_longer(
    cols = starts_with("n_"), # Select columns that start with "n_"
    names_to = "Timepoint", # The names of these columns will go into a new column "Timepoint"
    values_to = "Value" # The values will go into a new column "Value"
  ) %>%
  mutate(Timepoint = str_extract(Timepoint, "\\d+"), # Extract numeric part from Timepoint values
         Timepoint = as.numeric(Timepoint)) # Convert Timepoint values to numeric

df_long

# Plotting
ggplot(down_df_long, aes(x = Timepoint, y = Value, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = "ZM Treatment: downregulated genes",
       x = "time(hrs)",
       y = expression(paste(log[2](x), 'Change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
            plot.title.position = "plot")
  # geom_text(aes(label=symbol), hjust=0, vjust=0)
  # theme(legend.position = "none") # Remove legend to avoid overcrowding if you have many genes



down_df_long
