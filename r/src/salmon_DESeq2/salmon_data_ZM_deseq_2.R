getwd()
source("r/src/utils.R")
source("r/src/pca_utils.R")

treatment <- "ZM"
salmon_data_directory = file.path(getwd(), glue('data/organised/{treatment}/output_salmon'))
times = c(0, 16, 24, 36, 48)

print(salmon_data_directory)
# create_dds <- function(treatment, salmon_data_directory, times, file_prefix, replicates_list, batch_correction=TRUE, trim_data=TRUE) {
dds <-
  create_dds('ZM', salmon_data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "_data.RData")))


results <- results(dds)
resultsNames(dds)


results_t0_to_t16 <-
  results(dds, contrast = c("timepoint", "t16", "t0"))
results_t0_to_t16

results_t0_to_t48 <-
  results(
    dds,
    alpha = 0.05,
    lfcThreshold = 1,
    contrast = c("timepoint", "t48", "t0")
  )

results_t0_to_t48


summary(results_t0_to_t48)
plotMA(rest0t24)



# Other downstream analysis
# In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data.

normalised_data <- normTransform(dds)
library("vsn")
meanSdPlot(assay(normalised_data))

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced=~batch)
results_lrt <- results(dds_lrt)
results_lrt
summary(results_lrt)
# Shrinkage with ashr

resAshT <- lfcShrink(dds, coef=2, type="ashr", lfcThreshold=1)
plotMA(resAshT, ylim=c(-3,3), cex=.8)
# abline(h=c(-1,1), col="dodgerblue", lwd=2)



# Testing with thresholds

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()



## adding annotations
annotation <- get_annotation()
res = results(dds)
res = add_annotations_to_results(res)

annotated_results <- add_annotations_to_results(res)

resOrdered <- res[order(annotated_results$pvalue),]
resOrdered

summary(resOrdered)
mcols(res, use.names = TRUE)

res_05 = results(dds, alpha=0.05, contrast=c("timepoint", "t48", "t0"))
res_05
res_05 = add_annotations_to_results(res_05)
res_05
summary(res)
table(res$padj < 0.05)


# SEE http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#differential-expression-analysis
# How many genes were affected by treatement?
sum(res$pvalue < 0.05, na.rm=TRUE)


sum(res$padj < 0.1, na.rm=TRUE)


# Subset to determine the genes which have the strongest DOWN regulation.
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

# Subset to determine the genes which have the strongest UP regulation.
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])


topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("timepoint"))


hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")



library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("timepoint")])
pheatmap(mat, annotation_col = anno)
