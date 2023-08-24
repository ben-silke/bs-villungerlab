
library(EnhancedVolcano)
library(msigdbr)
setwd("~/bs-villungerlab/")
source("~/bs-villungerlab/r/src/pca_utils.R")
source("~/bs-villungerlab/r/src/utils.R")
source("~/bs-villungerlab/r/src/star_utils.R")
source("~/bs-villungerlab/r/src/gene_selection_utils.R")

source("~/bs-villungerlab/r/src/qPCR_utils.R")

zm_data <- read.csv('results/output_encode/ZM/all_ZM_gene_regulation_data.csv')

# zm_data <- read_csv("results/output_encode/ZM/fgsea_enrichment/fgsea_enrichment_ZM_increase_2n.csv")
View(zm_data)
# res <- results(dds)
# resultsNames(dds)
# res <- add_annotations_to_results(res)
# selected_genes <- as.character(resOrdered$symbol)
pathway_names <- as.character(zm_data$pathway)
pathway_names

plot <- EnhancedVolcano(zm_data,
                lab = pathway_names,
                x = 'NES',
                y = 'padj',
                xlim = c(-3,3),
                ylab = expression(paste('-Log'[10],' (adj P)')),
                xlab = 'NES',
                title = 'Pathway Enrichment ZM (2n)',
                pCutoff = 0.5,
                # at least double, or less than half
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)


ggsave(plot, filename = "results/output_encode/ZM/fgsea_enrichment/zm_2n_volcano_pathway_enrichment.pdf", width = 10, height = 10, limitsize = FALSE)

noc_data <- read_csv("results/output_encode/Noc/fgsea_enrichment/fgsea_enrichment_Noc_increase_2n.csv")
View(noc_data)
# res <- results(dds)
# resultsNames(dds)
# res <- add_annotations_to_results(res)
# selected_genes <- as.character(resOrdered$symbol)
noc_pathway_names <- as.character(noc_data$pathway)
noc_pathway_names

noc_plot <- EnhancedVolcano(noc_data,
                        lab = noc_pathway_names,
                        x = 'NES',
                        y = 'padj',
                        xlim = c(-3,3),
                        ylab = expression(paste('-Log'[10],' (adj P)')),
                        xlab = 'NES',
                        title = 'Pathway Enrichment Nocodazole (2n)',
                        pCutoff = 0.5,
                        # at least double, or less than half
                        FCcutoff = 1,
                        pointSize = 3.0,
                        labSize = 3.0)


ggsave(noc_plot, filename = "results/output_encode/Noc/fgsea_enrichment/noc_2n_volcano_pathway_enrichment.pdf", width = 10, height = 10, limitsize = FALSE)


