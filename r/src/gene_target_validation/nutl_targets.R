
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
library(plotly)
library(htmlwidgets)

setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

dir.create("results/output_encode/Nutl/target_genes")
dpi = 500
width_in <- 20  
height_in <- 20

load('~/bs-villungerlab/results/output_encode_1to6/Nutl_star_data.RData')
dds_Nutl <- ddseq_Nutl

data_file = "~/bs-villungerlab/results/output_encode_1to6/Nutl_unfiltered_results_files.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  results_Nutl_t8 <- get_unfiltered_results(dds_Nutl, "timepoint_t8_vs_t0", "_8")
  results_Nutl_t12 <- get_unfiltered_results(dds_Nutl, "timepoint_t12_vs_t0", "_12")
  results_Nutl_t16 <- get_unfiltered_results(dds_Nutl, "timepoint_t16_vs_t0", "_16")
  results_Nutl_t24 <- get_unfiltered_results(dds_Nutl, "timepoint_t24_vs_t0", "_24")
  results_Nutl_t48 <- get_unfiltered_results(dds_Nutl, "timepoint_t48_vs_t0", "_48")
  save(results_Nutl_t48, results_Nutl_t8, results_Nutl_t12, results_Nutl_t16, results_Nutl_t24, file=data_file)
}

colnames(results_Nutl_t8)
rownames(results_Nutl_t8)

all_df_merged_df <- merge_all_data(results_Nutl_t48, results_Nutl_t8, results_Nutl_t12, results_Nutl_t16, results_Nutl_t24, 'results/output_encode/Nutl/unfiltered_apeglm_Nutl_data.csv', 'full_join')
colnames(all_df_merged_df)
file_increase <- "results/output_encode/Nutl_full_signature_5n_increase_default_edge.csv"
file_decrease <- "results/output_encode/Nutl_full_signature_5n_decrease_default_edge.csv"

table_increase <- read.table(file_increase, sep=',', header=TRUE)
table_increase$Target.Gene
dim(table_increase)

merged_df <- all_df_merged_df

merged_df$symbol <- merged_df$symbol_48
merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_16[is.na(merged_df$symbol)]
merged_df$symbol[is.na(merged_df$symbol)] <- merged_df$symbol_8[is.na(merged_df$symbol)]

ls <- intersect(table_increase$Target.Gene, merged_df$symbol)

subset_df <- all_df_merged_df[merged_df$symbol %in% table_increase$Target.Gene, ]
subset_df$symbol_48
dim(subset_df)

subset_df_sorted_increase <- subset_df[order(-subset_df$log2FoldChange_16), ]
# set the number of results which you want
subset_df_sorted_increase_reduced <- head(subset_df_sorted_increase, 10)
subset_df_sorted_increase_reduced

# long_subset_df_sorted <- generate_complete_long_df(subset_df_sorted, 16)
# long_subset_df_sorted_plot <- plot_longdf(long_subset_df_sorted, "Nutl core up regulated genes | n10")
# long_subset_df_sorted_plot <- long_subset_df_sorted_plot + geom_text(aes(label = sprintf("p=%.3f", padj)), vjust = -1)
# long_subset_df_sorted_plot

# Functions
#######

df_long_increase_reduced <- generate_complete_long_df(subset_df_sorted_increase_reduced, 16)

# long_subset_df_sorted_plot <- long_subset_df_sorted_plot + geom_text(aes(label = sprintf("p=%.3f", padj)), vjust = -1)
plot_title <- "Nutl Core Genes: Increase | n10"
p_increase_reduced <- ggplot(df_long_increase_reduced, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y=expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "bottom")

p_increase_reduced
ggsave(filename = "results/output_encode/Nutl/target_genes/Nutl_p_increase_reduced.pdf", plot = p_increase_reduced, dpi=dpi, width=width_in, height=height_in)

### DECREASE
table_decrease <- read.table(file_decrease, sep=',', header=TRUE)
table_decrease$Target.Gene
dim(table_decrease)

merged_df_decrease <- all_df_merged_df

merged_df_decrease$symbol <- merged_df_decrease$symbol_48
merged_df_decrease$symbol[is.na(merged_df_decrease$symbol)] <- merged_df_decrease$symbol_16[is.na(merged_df_decrease$symbol)]
merged_df_decrease$symbol[is.na(merged_df_decrease$symbol)] <- merged_df_decrease$symbol_8[is.na(merged_df_decrease$symbol)]

ls <- intersect(table_decrease$Target.Gene, merged_df_decrease$symbol)

length(ls)

subset_df_decrease <- all_df_merged_df[merged_df_decrease$symbol %in% table_decrease$Target.Gene, ]
subset_df_decrease$symbol_48
dim(subset_df_decrease)

subset_df_decrease_sorted <- subset_df_decrease[order(subset_df_decrease$log2FoldChange_16), ]
# set the number of results which you want
subset_df_decrease_sorted_reduced <- head(subset_df_decrease_sorted, 10)
subset_df_decrease_sorted_reduced

#######
df_long_decrease <- generate_complete_long_df(subset_df_decrease_sorted_reduced, 16)

plot_title <- "Nutl Core Genes: n10"
p_decrease<- "he;"
p_decrease<- ggplot(df_long_decrease, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y=expression(paste(log[2](x), ' fold change'))) +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "bottom")

p_decrease
ggsave(filename = "results/output_encode/Nutl/target_genes/Nutl_p_decrease_reduced.pdf", plot = p_decrease, dpi=dpi, width=width_in, height=height_in)


# All Genes
######

df_long_decrease_all <- generate_complete_long_df(subset_df_decrease_sorted, 16)

plot_title <- "Nutl core genes: all downregulated"
p_decrease_all <- "he;"
p_decrease_all <- ggplot(df_long_decrease_all, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y='log 2 fold change') +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "bottom")

p_decrease_all
ggsave(filename = "results/output_encode/Nutl/target_genes/Nutl_p_decrease_all_reduced.pdf", plot = p_decrease_all, dpi=dpi, width=width_in, height=height_in)


df_long_increase <- generate_complete_long_df(subset_df_sorted_increase, 16)

plot_title <- "Nutl Core Genes"
p_increase_all<- "he;"
p_increase_all<- ggplot(df_long_increase, aes(x = Timepoint, y = log2foldchange, shape = symbol, group = symbol, color = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = plot_title,
       x = "time(hrs)",
       y='log 2 fold change') +
  theme(plot.title = element_text(hjust = 0.5), # Center the title
        plot.title.position = "plot",
        legend.position = "none")

p_increase_all
ggsave(filename = "results/output_encode/Nutl/target_genes/Nutl_p_increase_all.pdf", plot = p_increase_all, dpi=dpi, width=width_in, height=height_in)
 


##### INTERACTIVE PLOTS
##########
library(plotly)
library(htmlwidgets)

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
saveWidget(interactive_increase_plot, "results/output_encode/Nutl/target_genes/Nutl_interactive_increase_plot.html")

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
saveWidget(interactive_decrease_plot, "results/output_encode/Nutl/target_genes/Nutl_interactive_decrease_plot.html")



gene_data_file = "~/bs-villungerlab/results/output_encode_1to6/Nutl_fgsea_genes.RData"
load(gene_data_file)


subset_df_increase <- all_df_merged_df[merged_df$symbol %in% table_increase$Target.Gene, ]
subset_df_increase <- subset_df_increase[subset_df_increase$symbol %in% upregulated_genes, ]

View(subset_df_increase)

subset_df_decrease <- all_df_merged_df[merged_df_decrease$symbol %in% table_decrease$Target.Gene, ]


upregulated_genes_vec <- unlist(upregulated_genes)
downregulated_genes_vec <- unlist(downregulated_genes)
increase_shared <- intersect(table_increase$Target.Gene, upregulated_genes_vec)
increase_shared
decrease_shared <- intersect(table_decrease$Target.Gene, downregulated_genes_vec)
decrease_shared

subset_df_increase <- all_df_merged_df[merged_df$symbol %in% increase_shared, ]

View(subset_df_increase)
subset_df_increase$symbol_48

subset_df_decrease <- all_df_merged_df[merged_df$symbol %in% decrease_shared, ]

df_long_increase <- generate_complete_long_df(subset_df_increase, 16)
df_long_decrease <- generate_complete_long_df(subset_df_decrease, 16)

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
saveWidget(interactive_increase_plot, "results/output_encode/Nutl/merged/Nutl_interactive_increase_plot.html")


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
saveWidget(interactive_decrease_plot, "results/output_encode/Nutl/merged/Nutl_interactive_decrease_plot.html")


Nutl_df <- fix_labels(all_df_merged_df)


"FOXM1" %in% Nutl_df$symbol 
targets = c('BMF', "FOXM1", "SQSTM1", "NINJ1", "ZMAT3", "PHLDA3", "CCNA2", "CDCA8", "CDC25A", "AURKB", "ARID5B", "ARID5B", "ANKRD1")
subset_targets <- subset(Nutl_df, symbol %in% targets)
subset_targets

df_long_subset_targets <- generate_complete_long_df(subset_targets, 16)
df_long_subset_targets_plot <- plot_ly(
  df_long_subset_targets,
  x = ~Timepoint,
  y = ~log2foldchange,
  mode = "markers+lines",
  hoverinfo = "text",
  text = ~paste("T: ", Timepoint, "<br>l2fc: ", log2foldchange, "<br>padj: ", padj, '<br>ID: ', symbol),
  split= ~symbol
) %>%
  layout(title="Nutl Treatment: BMF + FOXM1")

df_long_subset_targets_plot
saveWidget(df_long_subset_targets_plot, "results/output_encode/Nutl/Nutl_bmf_foxm1.html")

