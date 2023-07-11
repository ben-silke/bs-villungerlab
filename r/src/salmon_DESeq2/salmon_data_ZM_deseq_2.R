library("vsn")
library("genefilter")


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



#######
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
######


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



#####
# Other downstream analysis
# In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data.

normalised_data <- normTransform(dds)
meanSdPlot(assay(normalised_data))

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
colData(vsd)

meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Demonstration of the variance increasing with mean
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)


##comparison between different normalisation methods
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

########



######
# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced=~batch)
results_lrt <- results(dds_lrt)
results_lrt
summary(results_lrt)
#######


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
annotation <- get_annotation(dds)
res = results(dds)
res = add_annotations_to_results(res)

annotated_results <- add_annotations_to_results(res)

# order results by padj value (most significant to least)
res <- res[order(res$padj),]

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



topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("timepoint")])
pheatmap(mat, annotation_col = anno)


######### Account for batch effects
library("RUVSeq")

set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)

pData(set)

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$batch, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}

### end account for batch effects