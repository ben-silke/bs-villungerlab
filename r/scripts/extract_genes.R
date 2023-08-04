
library(readxl)
EXTENSION = "lukas_proteomics/Omics/Proteomics/PD_exported/Reprocessed_3-pat1-vs-ctrl2&4_pat2-vs-ctrl1&3/"
# FILE_NAME = "CC_AVi_172-M1158-F1-F5-P13031-1-OTITOT-SumInd_TotalPeptideNorm-16plex-pat2_vs_control13-L2FC0.4_downreg_AR-20230418-V3.xlsx"
# SORTING_COLUMN = 'Abundance Ratio (log2): (DREAM-PL patient 2) / (Control groups 1 and 3)'
# P_VALUE_COL = 'Abundance Ratio Adj. P-Value: (DREAM-PL patient 2) / (Control groups 1 and 3)'


# FILE_NAME = "CC_AVi_172-M1158-F1-F5-P13031-1-OTITOT-SumInd_TotalPeptideNorm-16plex-pat1_vs_control24-L2FC0.4_downreg_AR-20230418-V3.xlsx"
FILE_NAME = "CC_AVi_172-M1158-F1-F5-P13031-1-OTITOT-SumInd_TotalPeptideNorm-16plex-pat1_vs_control24-L2FC0.4_upreg_AR-20230418-V3.xlsx"

SORTING_COLUMN = 'Abundance Ratio (log2): (DREAM-PL patient 1) / (Control groups 2 and 4)'
P_VALUE_COL = 'Abundance Ratio Adj. P-Value: (DREAM-PL patient 1) / (Control groups 2 and 4)'


GENE_SYMBOL = "Gene Symbol"

df = read_excel(file.path(EXTENSION, FILE_NAME))
df
dim(df)

names(df)[names(df) == SORTING_COLUMN] <- 'abundance_ratio_log2'
names(df)[names(df) == P_VALUE_COL] <- 'abundance_ratio_padj'
names(df)[names(df) == GENE_SYMBOL] <- 'gene_symbol'

colnames(df)


df <- subset(df, abundance_ratio_padj<0.05)
sorted_df <- df[order(-df$abundance_ratio_log2), ]  ##Upregulation
# sorted_df <- df[order(df$abundance_ratio_log2), ] ### down regulation

head(sorted_df)


# View(sorted_df)
dim(sorted_df)
sorted_df$abundance_ratio_log2
# sorted_df$abundance_ratio_padj

relevant_gene_list <- data.frame(symbol = sorted_df$gene_symbol)
colnames(relevant_gene_list)


noquote_gene_list <- noquote(relevant_gene_list$symbol)
noquote_gene_list <- c('gene_name', noquote_gene_list)

# write(noquote_gene_list, file = "results/lukas/p1c2c4_downr_relevant_genes_iregulon.txt")
write(noquote_gene_list, file = "results/lukas/p1c2c4_upr_relevant_genes_iregulon.txt")

