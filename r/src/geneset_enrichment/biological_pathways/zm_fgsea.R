#####
library(DESeq2)
library("apeglm")
library(org.Hs.eg.db)
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(htmlwidgets)
library("fgsea")
library(msigdbr)


setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")
######
load('~/bs-villungerlab/results/output_encode_1to6/ZM_star_data.RData')
dir.create("results/output_encode/ZM/fgsea_enrichment")


data_file = "~/bs-villungerlab/results/output_encode_1to6/ZM_unfiltered_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_ZM_t16 <- get_unfiltered_results(dds_ZM, "timepoint_t16_vs_t0", "_16")
  results_ZM_t20 <- get_unfiltered_results(dds_ZM, "timepoint_t20_vs_t0", "_20")
  results_ZM_t24 <- get_unfiltered_results(dds_ZM, "timepoint_t24_vs_t0", "_24")
  results_ZM_t36 <- get_unfiltered_results(dds_ZM, "timepoint_t36_vs_t0", "_36")
  results_ZM_t48 <- get_unfiltered_results(dds_ZM, "timepoint_t48_vs_t0", "_48")
  save(results_ZM_t48, results_ZM_t16, results_ZM_t20, results_ZM_t24, results_ZM_t36, file=data_file)
}

all_df_merged_df <- merge_all_data(results_ZM_t48, results_ZM_t16, results_ZM_t20, results_ZM_t24, results_ZM_t36, 'results/output_encode/ZM/all_ZM_gene_regulation_data.csv', 'full_join')
df <- fix_labels(all_df_merged_df)

data_file = "~/bs-villungerlab/results/output_encode_1to6/misgdrbr_df.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  msigdbr_df = msigdbr(species = "human", category="H")
  save(msigdbr_df, file=data_file)
}

msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

abs_foldchange_increase <- df[
  any(abs(df$log2FoldChange_16)>1 | 
  abs(df$log2FoldChange_20)>1 | 
  abs(df$log2FoldChange_24)>1 | 
  abs(df$log2FoldChange_36)>1 |
  abs(df$log2FoldChange_48)>1), ]


# include only 3 consecutive timepoints.
subset_df_3n <- abs_foldchange_increase[(!is.na(abs_foldchange_increase$log2FoldChange_16) & !is.na(abs_foldchange_increase$log2FoldChange_20) & !is.na(abs_foldchange_increase$log2FoldChange_24) |
  !is.na(abs_foldchange_increase$log2FoldChange_20) & !is.na(abs_foldchange_increase$log2FoldChange_24) & !is.na(abs_foldchange_increase$log2FoldChange_36) |
  !is.na(abs_foldchange_increase$log2FoldChange_24) & !is.na(abs_foldchange_increase$log2FoldChange_36) & !is.na(abs_foldchange_increase$log2FoldChange_48)), ]

# include only 2 consecutive timepoints.
subset_df_2n <- abs_foldchange_increase[
  !is.na(abs_foldchange_increase$log2FoldChange_16) & !is.na(abs_foldchange_increase$log2FoldChange_20) |
    !is.na(abs_foldchange_increase$log2FoldChange_20) & !is.na(abs_foldchange_increase$log2FoldChange_24) |
    !is.na(abs_foldchange_increase$log2FoldChange_24) & !is.na(abs_foldchange_increase$log2FoldChange_36) |
    !is.na(abs_foldchange_increase$log2FoldChange_36) & !is.na(abs_foldchange_increase$log2FoldChange_48)
  , ]

# Create Upper
######
# Need to sort the dataset before we provide to fgsea otherwise it will break.

## 3 
sorted_subset_foldchange_increase_3n <- subset_df_3n[order(-subset_df_3n$log2FoldChange_24), ]


sorted_foldchange_increase_3n_tovec <- data.frame(gene_symbol = sorted_subset_foldchange_increase_3n$symbol, log2foldchange <- sorted_subset_foldchange_increase_3n$log2FoldChange_24)
sorted_foldchange_increase_3n_tovec <- subset(sorted_foldchange_increase_3n_tovec, !is.na(log2foldchange))
nvec_increase_3n <- setNames(sorted_foldchange_increase_3n_tovec$log2foldchange, sorted_foldchange_increase_3n_tovec$gene_symbol)


fgseaRes_increase_3n <- fgsea(pathways = msigdbr_list,
                              stats    = nvec_increase_3n,
                              minSize  = 15,
                              maxSize  = 500)


fgseaRes_increase_3n_print <- fgseaRes_increase_3n
fgseaRes_increase_3n_print$leadingEdge <- sapply(fgseaRes_increase_3n_print$leadingEdge, toString)
write.csv(fgseaRes_increase_3n_print, file = "results/output_encode/ZM/fgsea_enrichment/fgsea_enrichment_ZM_increase_3n.csv")


## 2
sorted_subset_foldchange_increase_2n <- subset_df_2n[order(-subset_df_2n$log2FoldChange_24), ]

sorted_foldchange_increase_2n_tovec <- data.frame(gene_symbol = sorted_subset_foldchange_increase_2n$symbol, log2foldchange <- sorted_subset_foldchange_increase_2n$log2FoldChange_24)
sorted_foldchange_increase_2n_tovec <- subset(sorted_foldchange_increase_2n_tovec, !is.na(log2foldchange))
nvec_increase_2n <- setNames(sorted_foldchange_increase_2n_tovec$log2foldchange, sorted_foldchange_increase_2n_tovec$gene_symbol)

fgseaRes_increase_2n <- fgsea(pathways = msigdbr_list,
                              stats    = nvec_increase_2n,
                              minSize  = 15,
                              maxSize  = 500)

fgseaRes_increase_2n_print <- fgseaRes_increase_2n
fgseaRes_increase_2n_print$leadingEdge <- sapply(fgseaRes_increase_2n_print$leadingEdge, toString)
write.csv(fgseaRes_increase_2n_print, file = "results/output_encode/ZM/fgsea_enrichment/fgsea_enrichment_ZM_increase_2n.csv")


## ALSO - subset the pathways that have adj-p value <0.05, in which case we can use these to look at genes

# View(fgseaRes_increase_3n)
significant_pathways_3n <- subset(fgseaRes_increase_3n, 0.05>padj)

# head(significant_pathways, 3)
significant_pathways_2n <- subset(fgseaRes_increase_2n, 0.05>padj)
# View(significant_pathways_2n)

View(significant_pathways_3n)
significant_pathways_3n_print <- significant_pathways_3n
significant_pathways_3n_print$leadingEdge <- sapply(significant_pathways_3n_print$leadingEdge, toString)
write.csv(significant_pathways_3n_print, file = "results/output_encode/ZM/fgsea_enrichment/significant_pathways_3n.csv")

# sorted_subset_foldchange_increase_3n <- subset_df_3n[order(-subset_df_3n$log2FoldChange_24), ]

# Using 6 because we use top 3 motifs + top 3 TF's for the iRegulon Analysis
# significant_pathways_2n <- subset(fgseaRes_increase_2n, 0.05>padj)

upregulated_significant_pathways = subset(significant_pathways_3n, 0<NES)
downregulated_significant_pathways = subset(significant_pathways_3n, 0>NES)
View(top_significant_pathways)


upregulated_genes_list = upregulated_significant_pathways$leadingEdge
downregulated_genes_list = downregulated_significant_pathways$leadingEdge


upregulated_genes <- list()
for (gene_list in upregulated_genes_list) {
  for (gene in gene_list) {
    upregulated_genes <- append(upregulated_genes, gene)
  }
}

upregulated_genes

downregulated_genes <- list()
for (gene_list in downregulated_genes_list) {
  for (gene in gene_list) {
    downregulated_genes <- append(downregulated_genes, gene)
  }
}

gene_data_file = "~/bs-villungerlab/results/output_encode_1to6/ZM_fgsea_genes.RData"
save(upregulated_genes, downregulated_genes, file=gene_data_file)

# upregulated_genes
# downregulated_genes
downregulated_genes

subset_df_upregulated <- df[df$symbol %in% upregulated_genes, ]
subset_df_downregulated <- df[df$symbol %in% downregulated_genes, ]
View(subset_df_downregulated)

df_long_increase <- generate_complete_long_df(subset_df_upregulated, 24)
df_long_decrease <- generate_complete_long_df(subset_df_downregulated, 24)


interactive_increase_plot <- plot_ly(
  df_long_increase,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
)


interactive_increase_plot
saveWidget(interactive_increase_plot, "results/output_encode/ZM/fgsea_enrichment/ZM_interactive_increase_plot.html")


interactive_decrease_plot <- plot_ly(
  df_long_decrease,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
)

interactive_decrease_plot
saveWidget(interactive_decrease_plot, "results/output_encode/ZM/fgsea_enrichment/ZM_interactive_decrease_plot.html")



# fgseaRes_increase_2n <- fgsea(pathways = msigdbr_list,
#                               stats    = nvec_increase_2n,
#                               minSize  = 15,
#                               maxSize  = 500)

plotEnrichment(msigdbr_list[["HALLMARK_APOPTOSIS"]],
               nvec_increase_3n) + labs(title="Apoptosis")

plotEnrichment(msigdbr_list[["HALLMARK_P53_PATHWAY"]],
               nvec_increase_3n) + labs(title="p53")

# HALLMARK_G2M_CHECKPOINT
# HALLMARK_MYC_TARGETS_V1


plotEnrichment(msigdbr_list[["HALLMARK_G2M_CHECKPOINT"]],
               nvec_increase_3n) + labs(title="G2M Checkpoint")

plotEnrichment(msigdbr_list[["HALLMARK_MYC_TARGETS_V1"]],
               nvec_increase_3n) + labs(title="MYC Targets")
