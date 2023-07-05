## ---- echo=FALSE--------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE,message=FALSE)

## -----------------------------------------------------------------------------
library(tximportData)
dir <- system.file("extdata", package="tximportData")
list.files(dir)

## -----------------------------------------------------------------------------
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
all(file.exists(files))

## -----------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

## -----------------------------------------------------------------------------
library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)

## -----------------------------------------------------------------------------
library(tximport)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi$counts)

## -----------------------------------------------------------------------------
txi.tx <- tximport(files, type="salmon", txOut=TRUE)

## -----------------------------------------------------------------------------
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)

## -----------------------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene)
head(txi.salmon$counts)

## -----------------------------------------------------------------------------
tx2knownGene <- read_csv(file.path(dir, "tx2gene.csv"))
files <- file.path(dir,"sailfish", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
txi.sailfish <- tximport(files, type="sailfish", tx2gene=tx2knownGene)
head(txi.sailfish$counts)

## ----eval=FALSE---------------------------------------------------------------
#  txi <- tximport("quant.sf", type="none", txOut=TRUE,
#                  txIdCol="Name", abundanceCol="TPM",
#                  countsCol="NumReads", lengthCol="Length",
#                  importer=function(x) read_tsv(x, skip=8))

## -----------------------------------------------------------------------------
files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
txi.inf.rep <- tximport(files, type="salmon", txOut=TRUE)
names(txi.inf.rep)
names(txi.inf.rep$infReps)
dim(txi.inf.rep$infReps$sample1)

## -----------------------------------------------------------------------------
files <- file.path(dir, "kallisto_boot", samples$run, "abundance.h5")
names(files) <- paste0("sample",1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)
head(txi.kallisto$counts)

## -----------------------------------------------------------------------------
names(txi.kallisto)
names(txi.kallisto$infReps)
dim(txi.kallisto$infReps$sample1)

## -----------------------------------------------------------------------------
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
names(files) <- paste0("sample",1:6)
txi.kallisto.tsv <- tximport(files, type="kallisto", tx2gene=tx2gene,
                             ignoreAfterBar=TRUE)
head(txi.kallisto.tsv$counts)

## -----------------------------------------------------------------------------
files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results.gz"))
names(files) <- paste0("sample",1:6)
txi.rsem <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)
head(txi.rsem$counts)

## -----------------------------------------------------------------------------
files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".isoforms.results.gz"))
names(files) <- paste0("sample",1:6)
txi.rsem <- tximport(files, type="rsem", txIn=TRUE, txOut=TRUE)
head(txi.rsem$counts)

## ---- eval=FALSE--------------------------------------------------------------
#  tmp <- read_tsv(files[1])
#  tx2gene <- tmp[,c("t_name","gene_name")]
#  txi <- tximport(files, type="stringtie", tx2gene=tx2gene)

## ---- eval=FALSE--------------------------------------------------------------
#  files <- "path/to/alevin/quants_mat.gz"
#  txi <- tximport(files, type="alevin")

## ---- results="hide", messages=FALSE------------------------------------------
library(DESeq2)

## -----------------------------------------------------------------------------
sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=3)))
rownames(sampleTable) <- colnames(txi$counts)

## -----------------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
# dds is now ready for DESeq()
# see DESeq2 vignette

## ---- results="hide", messages=FALSE------------------------------------------
library(edgeR)

## -----------------------------------------------------------------------------
cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length,
# adjusted to avoid changing the magnitude of the counts.
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- cts / normMat

# Computing effective library sizes from scaled counts,
# to account for composition biases between samples.
library(edgeR)
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors,
# and calculating offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# filtering using the design information
design <- model.matrix(~ condition, data=sampleTable)
keep <- filterByExpr(y, design)
y <- y[keep,]
# y is now ready for estimate dispersion functions
# see edgeR User's Guide

## -----------------------------------------------------------------------------
cpms <- edgeR::cpm(y, offset=y$offset, log=FALSE)

## -----------------------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
txi <- tximport(files, type="salmon",
                tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)

# filtering using the design information:
design <- model.matrix(~ condition, data=sampleTable)
keep <- filterByExpr(y, design)
y <- y[keep,]

# normalize and run voom transformation
y <- calcNormFactors(y)
v <- voom(y, design)
# v is now ready for lmFit()
# see limma User's Guide

## -----------------------------------------------------------------------------
sessionInfo()
