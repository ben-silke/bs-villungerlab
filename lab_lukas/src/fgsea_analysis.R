
library(readxl)
library(plotly)
library(htmlwidgets)
library(dplyr)
library("fgsea")
library(EnhancedVolcano)

setwd("~/bs-villungerlab/")

data_file = "~/bs-villungerlab/references/misgdrbr_df.Rdata"
if (file.exists(data_file)){
  load(data_file)
  print("file_exists")
} else {
  msigdbr_df = msigdbr(species = "human", category="H")
  save(msigdbr_df, file=data_file)
}

msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

df <- read_excel("lab_lukas/Omics/RNAseq/DESeq2_1.18.1.DREAM_PL_VS_control.alpha05.xlsx", )
# View(df)
df <- df[order(-df$log2FoldChange.shrunk), ]

df_tovec <- data.frame(gene_symbol = df$gene_id, log2foldchange <- df$log2FoldChange.shrunk)
df_tovec <- subset(df_tovec, !is.na(log2foldchange))
df_tovec
vec <- setNames(df_tovec$log2foldchange, df_tovec$gene_symbol)
vec
dreampl_fgseaRes <- fgsea(pathways = msigdbr_list,
                              stats    = vec,
                              minSize  = 15,
                              maxSize  = 500)

dreampl_fgseaRes_print <- dreampl_fgseaRes
dreampl_fgseaRes_print$leadingEdge <- sapply(dreampl_fgseaRes_print$leadingEdge, toString)

# View(dreampl_fgseaRes)

write.csv(dreampl_fgseaRes_print, file="results/lab_lukas/geneset_enrichment.xlsx")


dreampl_pathwaynames <- as.character(dreampl_fgseaRes$pathway)
dreampl_pathwaynames

fgsea_plot <- EnhancedVolcano(dreampl_fgseaRes,
                        lab = dreampl_pathwaynames,
                        x = 'NES',
                        y = 'padj',
                        xlim = c(-3,3),
                        ylab = expression(paste('-Log'[10],' (adj P)')),
                        xlab = 'NES',
                        title = 'Dream PL pathway Enrichment',
                        pCutoff = 0.5,
                        # at least double, or less than half
                        FCcutoff = 1,
                        pointSize = 3.0,
                        labSize = 3.0)

fgsea_plot
ggsave(fgsea_plot, filename = "results/lab_lukas/fgsea_enrichment.pdf", width = 10, height = 10, limitsize = FALSE)
