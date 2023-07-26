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

# load('~/bs-villungerlab/results/output_encode_1to6/Noc_star_data.RData')
# dds_noc <- ddseq_Noc
# results_noc_t48Sig <- subset(results_noc_t48, padj < 0.1)
# results_noc_t48 <- lfcShrink(dds_noc, coef="timepoint_t48_vs_t0", type="apeglm")


results_zm_t16 <- lfcShrink(dds_zm, coef="timepoint_t16_vs_t0", type="apeglm")
results_zm_t16Sig <- subset(results_zm_t16, padj < 0.1)
decreasing_t16_zm <- head(results_zm_t16Sig[ order(results_zm_t16Sig$log2FoldChange), ], 10)
increasing_t16_zm <-head(results_zm_t16Sig[ order(results_zm_t16Sig$log2FoldChange, decreasing = TRUE), ], 10)
decreasing_t16_zm$timepoint <- 16
increasing_t16_zm$timepoint <- 16

results_zm_t20 <- lfcShrink(dds_zm, coef="timepoint_t20_vs_t0", type="apeglm")
results_zm_t20Sig <- subset(results_zm_t20, padj < 0.1)
decreasing_t20_zm <- head(results_zm_t20Sig[ order(results_zm_t20Sig$log2FoldChange), ], 10)
increasing_t20_zm <-head(results_zm_t20Sig[ order(results_zm_t20Sig$log2FoldChange, decreasing = TRUE), ], 10)
decreasing_t20_zm$timepoint <- 20
increasing_t20_zm$timepoint <- 20


results_zm_t24 <- lfcShrink(dds_zm, coef="timepoint_t24_vs_t0", type="apeglm")
results_zm_t24Sig <- subset(results_zm_t24, padj < 0.1)
decreasing_t24_zm <- head(results_zm_t24Sig[ order(results_zm_t24Sig$log2FoldChange), ], 10)
increasing_t24_zm <-head(results_zm_t24Sig[ order(results_zm_t24Sig$log2FoldChange, decreasing = TRUE), ], 10)
decreasing_t24_zm$timepoint <- 24
increasing_t24_zm$timepoint <- 24

results_zm_t36 <- lfcShrink(dds_zm, coef="timepoint_t36_vs_t0", type="apeglm")
results_zm_t36Sig <- subset(results_zm_t36, padj < 0.1)
decreasing_t36_zm <- head(results_zm_t36Sig[ order(results_zm_t36Sig$log2FoldChange), ], 10)
increasing_t36_zm <-head(results_zm_t36Sig[ order(results_zm_t36Sig$log2FoldChange, decreasing = TRUE), ], 10)
decreasing_t36_zm$timepoint <- 36
increasing_t36_zm$timepoint <- 36


results_zm_t48 <- lfcShrink(dds_zm, coef="timepoint_t48_vs_t0", type="apeglm")
results_zm_t48Sig <- subset(results_zm_t48, padj < 0.1)
decreasing_t48_zm <- head(results_zm_t48Sig[ order(results_zm_t48Sig$log2FoldChange), ], 10)
increasing_t48_zm <-head(results_zm_t48Sig[ order(results_zm_t48Sig$log2FoldChange, decreasing = TRUE), ], 10)
decreasing_t48_zm$timepoint <- 48
increasing_t48_zm$timepoint <- 48

colnames(increasing_t16_zm)
colnames(increasing_t20_zm)
colnames(increasing_t24_zm)
colnames(increasing_t36_zm)
colnames(increasing_t48_zm)

increasing_t16_zm <- add_annotations_to_results(increasing_t16_zm)
increasing_t20_zm <- add_annotations_to_results(increasing_t20_zm)
increasing_t24_zm <- add_annotations_to_results(increasing_t24_zm)
increasing_t36_zm <- add_annotations_to_results(increasing_t36_zm)
increasing_t48_zm <- add_annotations_to_results(increasing_t48_zm)

decreasing_t16_zm <- add_annotations_to_results(decreasing_t16_zm)
decreasing_t20_zm <- add_annotations_to_results(decreasing_t20_zm)
decreasing_t24_zm <- add_annotations_to_results(decreasing_t24_zm)
decreasing_t36_zm <- add_annotations_to_results(decreasing_t36_zm)
decreasing_t48_zm <- add_annotations_to_results(decreasing_t48_zm)


increasing_df <- rbind(increasing_t16_zm, increasing_t20_zm, increasing_t24_zm, increasing_t36_zm, increasing_t48_zm)

decreasing_df <- rbind(decreasing_t16_zm, decreasing_t20_zm, decreasing_t24_zm, decreasing_t36_zm, decreasing_t48_zm)
increasing_df <- as.data.frame(increasing_df)

# long_inc_df <- increasing_df %>% gather(key = "timepoint", value = "log2FoldChange", -symbol)
# long_inc_df <- decreasing_df %>% gather(key = "timepoint", value = "log2FoldChange", -symbol)


increasing_df <- as.data.frame(increasing_df)

ggplot(data = increasing_df, aes(x = timepoint, y = log2FoldChange, group = symbol, colour = symbol)) +
  geom_line() +
  labs(title = "Change over time", x = "Time Point", y = "Log2 Fold Change")


  