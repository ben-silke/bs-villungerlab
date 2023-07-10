#setting the environment
setwd("~/OneDrive - CeMM Research Center GmbH/Year 3/RNA sequencing/New Sequencing/p53_core_ZM_interactive/")
library(ggplot2)
library(ggrepel)
library(plotly)
library(htmlwidgets)
library(readxl)

#importing files from the facility
t16_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/hallmark_apoptosis_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t16_vs_ZM_t00_against_intercept_genes.tsv")
t20_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/hallmark_apoptosis_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t20_vs_ZM_t00_against_intercept_genes.tsv")
t24_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/hallmark_apoptosis_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t24_vs_ZM_t00_against_intercept_genes.tsv")
t36_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/hallmark_apoptosis_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t36_vs_ZM_t00_against_intercept_genes.tsv")
t48_all <- read.delim("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/RNA sequencing/New Sequencing/hallmark_apoptosis_ZM/rnaseq_deseq_ZM2_contrast_group_ZM_t48_vs_ZM_t00_against_intercept_genes.tsv")

#importing genes defined previously from fisher 2017 & andrysik 
core_p53 <- read_excel("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/TimeCourse Analysis/Background_Datasets/core_p53_andrysik.xlsx")
up_genes <- core_p53$`Gene name`
down_reg <- read.table("~/Library/CloudStorage/OneDrive-CeMMResearchCenterGmbH/Year 3/TimeCourse Analysis/Background_Datasets/p53-DREAM_down.txt", quote="\"", comment.char="")
down_genes <- down_reg$V1


###PLOTS###

##t16##
t16_noNA <- t16_all[complete.cases(t16_all$padj),]
t16_significant <- subset(t16_noNA, t16_noNA$padj <= 0.1)
t16_non_significant <- subset(t16_noNA, t16_noNA$padj > 0.1)
t16_significant$pathway <- ifelse(t16_significant$gene_name %in% up_genes, "p53-core", 
                               ifelse(t16_significant$gene_name %in% down_genes, "p53-DREAM_downreg", "no"))

t16_significant$regulation <- ifelse(t16_significant$log2FoldChange <= -1 | t16_significant$log2FoldChange >= 1, "Regulated", "Not regulated")
# Significant Plot
t16_significant_plot <- ggplot(t16_significant, aes(x = log2FoldChange, y = -log10(padj), text = gene_name, shape = regulation, color = pathway)) +
  geom_point(aes(color = pathway)) +
  scale_shape_manual(values = c("Regulated" = 16, "Not regulated" = 15)) +
  scale_color_manual(values = c("p53-core" = "#005564", "p53-DREAM_downreg" = "#F8B100", "no" = "grey"), name = "p53 responsive genes") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significant genes") +
  ggtitle("treatment: ZM - timepoint: 16h") +
  xlim(-5, 5) +
  guides(shape = guide_legend(title = "Regulation", override.aes = list(size = 3), labels = c("A","B")))

#t16_significant_plot

t16_significant_plot_interactive <- ggplotly(t16_significant_plot, tooltip = "text")

#t16_significant_plot_interactive

saveWidget(t16_significant_plot_interactive, "ZM_t16_significant.html", selfcontained = TRUE)


##t20##
t20_noNA <- t20_all[complete.cases(t20_all$padj),]
t20_significant <- subset(t20_noNA, t20_noNA$padj <= 0.1)
t20_non_significant <- subset(t20_noNA, t20_noNA$padj > 0.1)
t20_significant$pathway <- ifelse(t20_significant$gene_name %in% up_genes, "p53-core", 
                               ifelse(t20_significant$gene_name %in% down_genes, "p53-DREAM_downreg", "no"))
t20_significant$regulation <- ifelse(t20_significant$log2FoldChange <= -1 | t20_significant$log2FoldChange >= 1, "Regulated", "Not regulated")
# Sigificant Plot
t20_significant_plot <- ggplot(t20_significant, aes(x = log2FoldChange, y = -log10(padj), text = gene_name, shape = t20_significant$regulation, color = pathway)) +
  geom_point(aes(color = pathway)) +
  scale_shape_manual(values = c("Regulated" = 20, "Not regulated" = 15)) +
  scale_color_manual(values = c("p53-core" = "#005564", "p53-DREAM_downreg" = "#F8B100", "no" = "grey"), name = "p53 responsive genes") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significant genes") +
  ggtitle("treatment: ZM - timepoint: 20h") +
  xlim(-5, 5) +
  guides(shape = guide_legend(title = "Regulation", override.aes = list(size = 3), labels = c("A","B")))

#t20_significant_plot

t20_significant_plot_interactive <- ggplotly(t20_significant_plot, tooltip = "text")

saveWidget(t20_significant_plot_interactive, "ZM_t20_significant.html", selfcontained = TRUE)


##t24##
t24_noNA <- t24_all[complete.cases(t24_all$padj),]
t24_significant <- subset(t24_noNA, t24_noNA$padj <= 0.1)
t24_non_significant <- subset(t24_noNA, t24_noNA$padj > 0.1)
t24_significant$pathway <- ifelse(t24_significant$gene_name %in% up_genes, "p53-core", 
                                  ifelse(t24_significant$gene_name %in% down_genes, "p53-DREAM_downreg", "no"))
t24_significant$regulation <- ifelse(t24_significant$log2FoldChange <= -1 | t24_significant$log2FoldChange >= 1, "Regulated", "Not regulated")
# Sigificant Plot
t24_significant_plot <- ggplot(t24_significant, aes(x = log2FoldChange, y = -log10(padj), text = gene_name, shape = t24_significant$regulation, color = pathway)) +
  geom_point(aes(color = pathway)) +
  scale_shape_manual(values = c("Regulated" = 20, "Not regulated" = 15)) +
  scale_color_manual(values = c("p53-core" = "#005564", "p53-DREAM_downreg" = "#F8B100", "no" = "grey"), name = "p53 responsive genes") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significant genes") +
  ggtitle("treatment: ZM - timepoint: 24h") +
  xlim(-5, 5) +
  guides(shape = guide_legend(title = "Regulation", override.aes = list(size = 3), labels = c("A","B")))

t24_significant_plot

t24_significant_plot_interactive <- ggplotly(t24_significant_plot, tooltip = "text")

saveWidget(t24_significant_plot_interactive, "ZM_t24_significant.html", selfcontained = TRUE)



##t36##
t36_noNA <- t36_all[complete.cases(t36_all$padj),]
t36_significant <- subset(t36_noNA, t36_noNA$padj <= 0.1)
t36_non_significant <- subset(t36_noNA, t36_noNA$padj > 0.1)
t36_significant$pathway <- ifelse(t36_significant$gene_name %in% up_genes, "p53-core", 
                                  ifelse(t36_significant$gene_name %in% down_genes, "p53-DREAM_downreg", "no"))
t36_significant$regulation <- ifelse(t36_significant$log2FoldChange <= -1 | t36_significant$log2FoldChange >= 1, "Regulated", "Not regulated")
# Sigificant Plot
t36_significant_plot <- ggplot(t36_significant, aes(x = log2FoldChange, y = -log10(padj), text = gene_name, shape = t36_significant$regulation, color = pathway)) +
  geom_point(aes(color = pathway)) +
  scale_shape_manual(values = c("Regulated" = 20, "Not regulated" = 15)) +
  scale_color_manual(values = c("p53-core" = "#005564", "p53-DREAM_downreg" = "#F8B100", "no" = "grey"), name = "p53 responsive genes") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significant genes") +
  ggtitle("treatment: ZM - timepoint: 36h") +
  xlim(-5, 5) +
  guides(shape = guide_legend(title = "Regulation", override.aes = list(size = 3), labels = c("A","B")))

t36_significant_plot

t36_significant_plot_interactive <- ggplotly(t36_significant_plot, tooltip = "text")

saveWidget(t36_significant_plot_interactive, "ZM_t36_significant.html", selfcontained = TRUE)


##t48##
t48_noNA <- t48_all[complete.cases(t48_all$padj),]
t48_significant <- subset(t48_noNA, t48_noNA$padj <= 0.1)
t48_non_significant <- subset(t48_noNA, t48_noNA$padj > 0.1)
t48_significant$pathway <- ifelse(t48_significant$gene_name %in% up_genes, "p53-core", 
                                  ifelse(t48_significant$gene_name %in% down_genes, "p53-DREAM_downreg", "no"))
t48_significant$regulation <- ifelse(t48_significant$log2FoldChange <= -1 | t48_significant$log2FoldChange >= 1, "Regulated", "Not regulated")
# Sigificant Plot
t48_significant_plot <- ggplot(t48_significant, aes(x = log2FoldChange, y = -log10(padj), text = gene_name, shape = t48_significant$regulation, color = pathway)) +
  geom_point(aes(color = pathway)) +
  scale_shape_manual(values = c("Regulated" = 20, "Not regulated" = 15)) +
  scale_color_manual(values = c("p53-core" = "#005564", "p53-DREAM_downreg" = "#F8B100", "no" = "grey"), name = "p53 responsive genes") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  #theme(legend.position = "none") +
  labs(x = "Log2 Fold Change", y = "-log10(adj P-value)", color = "Significant genes") +
  ggtitle("treatment: ZM - timepoint: 48h") +
  xlim(-5, 5) +
  guides(shape = guide_legend(title = "Regulation", override.aes = list(size = 3), labels = c("A","B")))

t48_significant_plot

t48_significant_plot_interactive <- ggplotly(t48_significant_plot, tooltip = "text")

saveWidget(t48_significant_plot_interactive, "ZM_t48_significant.html", selfcontained = TRUE)
