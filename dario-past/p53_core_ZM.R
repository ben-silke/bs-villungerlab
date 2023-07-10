setwd("~/OneDrive - CeMM Research Center GmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM")

#importing files from the facility
t16_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t16_vs_ZM_t00_against_intercept_genes.tsv")
t20_all <-read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t20_vs_ZM_t00_against_intercept_genes.tsv")
t24_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t24_vs_ZM_t00_against_intercept_genes.tsv")
t36_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t36_vs_ZM_t00_against_intercept_genes.tsv")
t48_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t48_vs_ZM_t00_against_intercept_genes.tsv")

#importing the p53 core genes defined previously
library(readxl)
core_p53 <- read_excel("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/TimeCourse Analysis/Background_Datasets/core_p53_andrysik.xlsx")
up_genes <- core_p53$`Gene name`
down_reg <- read.table("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/TimeCourse Analysis/Background_Datasets/p53-DREAM_down.txt", quote="\"", comment.char="")
down_genes <- down_reg$V1

##t16##
library(ggplot2)
library(ggrepel)
t16_significant <- subset(t16_all, t16_all$padj <= 0.1)
t16_p53core_significant_labels <- subset(t16_significant, t16_significant$gene_name %in% core_p53$`Gene name`)
t16_p53signature_significant_labels <- subset(t16_significant, t16_significant$gene_name %in% core_p53$`Gene name` | t16_significant$gene_name %in% down_reg$V1)
#calculate the number of p53 targets over the upregulated DEGs
t16_number_up <- paste0((nrow(subset(t16_significant, t16_significant$gene_name %in% core_p53$`Gene name` 
                                     & t16_significant$log2FoldChange >= 1))), "/",
                        (nrow(subset(t16_significant, t16_significant$log2FoldChange >= 1))))
t16_number_down <- paste0((nrow(subset(t16_significant, t16_significant$gene_name %in% down_reg$V1 
                                       & t16_significant$log2FoldChange <= -1))), "/",
                          (nrow(subset(t16_significant, t16_significant$log2FoldChange <= -1))))

t16_vulcano_plot <- ggplot(t16_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= 1 & (-log10(padj)) >= 1, 
                                "significant", "not_significant")), alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(significant = "red", not_significant = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.title.align = 0.5) +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significance") +
  theme(axis.title.x = element_text(size = 9),  
        axis.title.y = element_text(size = 9)) +
  ggtitle("treatment: ZM - timepoint: 16h") +
  theme(plot.title = element_text(face = "bold", color = "blue", size = 16, hjust = 0.5)) +
  xlim(-6,6) +
  geom_text_repel(data = t16_p53signature_significant_labels, 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_name), 
                  box.padding = 0.1, point.padding = 0.5, force = 10, 
                  nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 15) + 
  annotate("text", x = 5 , y = 7.5, label = t16_number_up, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = 30)) +
  annotate("text", x = -5 , y = 7.5, label = t16_number_down, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = -30))

t16_vulcano_plot
ggsave("Output/t16_vulcano_plot.png", plot = t16_vulcano_plot, width = 8, height = 6, dpi = 300)


##t20##
t20_significant <- subset(t20_all, t20_all$padj <= 0.1)
#t20_p53core_significant_labels <- subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name`)
t20_p53signature_significant_labels <- subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name` | t20_significant$gene_name %in% down_reg$V1)
#calculate the number of p53 targets over the upregulated DEGs
t20_number_up <- paste0((nrow(subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name` 
                                     & t20_significant$log2FoldChange >= 1))), "/",
                        (nrow(subset(t20_significant, t20_significant$log2FoldChange >= 1))))
t20_number_down <- paste0((nrow(subset(t20_significant, t20_significant$gene_name %in% down_reg$V1 
                                       & t20_significant$log2FoldChange <= -1))), "/",
                          (nrow(subset(t20_significant, t20_significant$log2FoldChange <= -1))))

t20_vulcano_plot <- ggplot(t20_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= 1 & (-log10(padj)) >= 1, 
                                "significant", "not_significant")), alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(significant = "red", not_significant = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.title.align = 0.5) +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significance") +
  theme(axis.title.x = element_text(size = 9),  
        axis.title.y = element_text(size = 9)) +
  ggtitle("treatment: ZM - timepoint: 20h") +
  theme(plot.title = element_text(face = "bold", color = "blue", size = 16, hjust = 0.5)) +
  xlim(-6,6) +
  geom_text_repel(data = t20_p53signature_significant_labels, 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_name), 
                  box.padding = 0.1, point.padding = 0.5, force = 10, 
                  nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 15) + 
  annotate("text", x = 5 , y = 12.5, label = t20_number_up, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = 30)) +
  annotate("text", x = -5 , y = 12.5, label = t20_number_down, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = -30))

t20_vulcano_plot
ggsave("Output/t20_vulcano_plot.png", plot = t20_vulcano_plot, width = 8, height = 6, dpi = 300)


##t24##
t24_significant <- subset(t24_all, t24_all$padj <= 0.1)
#t20_p53core_significant_labels <- subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name`)
t24_p53signature_significant_labels <- subset(t24_significant, t24_significant$gene_name %in% core_p53$`Gene name` 
                                              | t24_significant$gene_name %in% down_reg$V1)
#calculate the number of p53 targets over the upregulated DEGs
t24_number_up <- paste0((nrow(subset(t24_significant, t24_significant$gene_name %in% core_p53$`Gene name` 
                                     & t24_significant$log2FoldChange >= 1))), "/",
                        (nrow(subset(t24_significant, t24_significant$log2FoldChange >= 1))))
t24_number_down <- paste0((nrow(subset(t24_significant, t24_significant$gene_name %in% down_reg$V1 
                                       & t24_significant$log2FoldChange <= -1))), "/",
                          (nrow(subset(t24_significant, t24_significant$log2FoldChange <= -1))))

t24_vulcano_plot <- ggplot(t24_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= 1 & (-log10(padj)) >= 1, 
                                "significant", "not_significant")), alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(significant = "red", not_significant = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.title.align = 0.5) +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significance") +
  theme(axis.title.x = element_text(size = 9),  
        axis.title.y = element_text(size = 9)) +
  ggtitle("treatment: ZM - timepoint: 24h") +
  theme(plot.title = element_text(face = "bold", color = "blue", size = 16, hjust = 0.5)) +
  xlim(-6,6) +
  geom_text_repel(data = t24_p53signature_significant_labels, 
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_name), 
                  box.padding = 0.1, point.padding = 0.5, force = 10, 
                  nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 15) + 
  annotate("text", x = 5 , y = 22, label = t24_number_up, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = 30)) +
  annotate("text", x = -5 , y = 22, label = t24_number_down, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = -30))
  
t24_vulcano_plot
ggsave("Output/t24vulcano_plot.png", plot = t24_vulcano_plot, width = 8, height = 6, dpi = 300)


##t36##
t36_significant <- subset(t36_all, t36_all$padj <= 0.1)
#t36_p53core_significant_labels <- subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name`)
t36_p53signature_significant_labels <- subset(t36_significant, t36_significant$gene_name %in% core_p53$`Gene name` | t36_significant$gene_name %in% down_reg$V1)
#calculate the number of p53 targets over the upregulated DEGs
t36_number_up <- paste0((nrow(subset(t36_significant, t36_significant$gene_name %in% core_p53$`Gene name` 
                                     & t36_significant$log2FoldChange >= 1))), "/",
                        (nrow(subset(t36_significant, t36_significant$log2FoldChange >= 1))))
t36_number_down <- paste0((nrow(subset(t36_significant, t36_significant$gene_name %in% down_reg$V1 
                                       & t36_significant$log2FoldChange <= -1))), "/",
                          (nrow(subset(t36_significant, t36_significant$log2FoldChange <= -1))))

t36_vulcano_plot <- ggplot(t36_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= 1 & (-log10(padj)) >= 1, "significant", "not_significant")), alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(significant = "red", not_significant = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",legend.title.align = 0.5) +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significance") +
  theme(axis.title.x = element_text(size = 9),  
        axis.title.y = element_text(size = 9)) +
  ggtitle("treatment: ZM - timepoint: 36h") +
  theme(plot.title = element_text(face = "bold", color = "blue", size = 16, hjust = 0.5)) +
  xlim(-6,6) +
  geom_text_repel(data = t36_p53signature_significant_labels, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), box.padding = 0.1, point.padding = 0.5, force = 10, nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 15) + 
  annotate("text", x = 5 , y = 35, label = t36_number_up, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = 30)) +
  annotate("text", x = -5 , y = 35, label = t36_number_down, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = -30))

t36_vulcano_plot
ggsave("Output/t36vulcano_plot.png", plot = t36_vulcano_plot, width = 8, height = 6, dpi = 300)


##t48##
t48_significant <- subset(t48_all, t48_all$padj <= 0.1)
#t36_p53core_significant_labels <- subset(t20_significant, t20_significant$gene_name %in% core_p53$`Gene name`)
t48_p53signature_significant_labels <- subset(t48_significant, t48_significant$gene_name %in% core_p53$`Gene name` | t48_significant$gene_name %in% down_reg$V1)
#calculate the number of p53 targets over the upregulated DEGs
t48_number_up <- paste0((nrow(subset(t48_significant, t48_significant$gene_name %in% core_p53$`Gene name` 
                                     & t48_significant$log2FoldChange >= 1))), "/",
                        (nrow(subset(t48_significant, t48_significant$log2FoldChange >= 1))))
t48_number_down <- paste0((nrow(subset(t48_significant, t48_significant$gene_name %in% down_reg$V1 
                                       & t48_significant$log2FoldChange <= -1))), "/",
                          (nrow(subset(t48_significant, t48_significant$log2FoldChange <= -1))))

t48_vulcano_plot <- ggplot(t48_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) >= 1 & (-log10(padj)) >= 1, "significant", "not_significant")), alpha = 0.8, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(significant = "red", not_significant = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",legend.title.align = 0.5) +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significance") +
  theme(axis.title.x = element_text(size = 9),  
        axis.title.y = element_text(size = 9)) +
  ggtitle("treatment: ZM - timepoint: 48h") +
  theme(plot.title = element_text(face = "bold", color = "blue", size = 16, hjust = 0.5)) +
  xlim(-6,6) +
  geom_text_repel(data = t48_p53signature_significant_labels, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), box.padding = 0.1, point.padding = 0.5, force = 10, nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 15) + 
  annotate("text", x = 5 , y = 55, label = t48_number_up, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = 30)) +
  annotate("text", x = -5 , y = 55, label = t48_number_down, hjust = 0, vjust = 0) + 
  theme(plot.margin = margin(b = -30))

t48_vulcano_plot

ggsave("Output/t48vulcano_plot.png", plot = t48_vulcano_plot, width = 8, height = 6, dpi = 300)



# > sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Ventura 13.4
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.3 ggplot2_3.4.2 readxl_1.4.2 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10       rstudioapi_0.14   magrittr_2.0.3    tidyselect_1.2.0  munsell_0.5.0    
# [6] colorspace_2.1-0  R6_2.5.1          ragg_1.2.5        rlang_1.1.1       fansi_1.0.4      
# [11] dplyr_1.1.2       tools_4.2.2       grid_4.2.2        gtable_0.3.3      utf8_1.2.3       
# [16] cli_3.6.1         withr_2.5.0       systemfonts_1.0.4 tibble_3.2.1      lifecycle_1.0.3  
# [21] textshaping_0.3.6 farver_2.1.1      vctrs_0.6.2       glue_1.6.2        labeling_0.4.2   
# [26] compiler_4.2.2    pillar_1.9.0      cellranger_1.1.0  generics_0.1.3    scales_1.2.1     
# [31] pkgconfig_2.0.3  
