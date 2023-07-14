library(EnhancedVolcano)

source("r/src/utils.R")
source("r/src/pca_utils.R")


getwd()
load("r/data/ZM_data.RData")


res <- results(dds)
resultsNames(dds)
res <- add_annotations_to_results(res)
selected_genes <- as.character(resOrdered$symbol)


EnhancedVolcano(resOrdered,
                lab = selected_genes,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8,8),
                title = 'Differential expression',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)


sum(resOrdered$padj < 0.05, na.rm=TRUE)
