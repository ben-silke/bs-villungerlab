library("airway")
library("DESeq2")



data("airway")
se_airway <- airway

dds <- DESeqDataSet(se_airway, design = ~ cell + dex)


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition

dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

# dds$condition <- relevel(dds$condition, ref = "untreated")



dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, name="condition_treated_vs_untreated")

# 
# res <- results(dds, contrast=c("condition","treated","untreated"))

resultsNames(dds)


# resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC <- lfcShrink(dds, coef="dex_untrt_vs_trt", type="apeglm")
resLFC


resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)


plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]