#tximeta install
#BiocManager::install("tximeta")
#explore the vignette
#browseVignettes("tximeta")

homedir <- ("/Users/dario/OneDrive - CeMM Research Center GmbH/Year 3/RNA sequencing/New Sequencing/Nutl_analysis/")
setwd(homedir)
quantdir <- file.path(homedir,"quantification")
resdir <- file.path(homedir, "Nutl_DESeq2")

#importing the different quant.sf files as indipendent variables
t0_r1 <- file.path(quantdir,"Nutl_t0_r1","quant.sf")
t0_r2 <- file.path(quantdir,"Nutl_t0_r2","quant.sf")
t0_r3 <- file.path(quantdir,"Nutl_t0_r3","quant.sf")

t8_r1 <- file.path(quantdir,"Nutl_t8_r1","quant.sf")
t8_r2 <- file.path(quantdir,"Nutl_t8_r2","quant.sf")
t8_r3 <- file.path(quantdir,"Nutl_t8_r3","quant.sf")

t12_r1 <- file.path(quantdir,"Nutl_t12_r1","quant.sf")
t12_r2 <- file.path(quantdir,"Nutl_t12_r2","quant.sf")
t12_r3 <- file.path(quantdir,"Nutl_t12_r3","quant.sf")

t16_r1 <- file.path(quantdir,"Nutl_t16_r1","quant.sf")
t16_r2 <- file.path(quantdir,"Nutl_t16_r2","quant.sf")
t16_r3 <- file.path(quantdir,"Nutl_t16_r3","quant.sf")

t24_r1 <- file.path(quantdir,"Nutl_t24_r1","quant.sf")
t24_r2 <- file.path(quantdir,"Nutl_t24_r2","quant.sf")
t24_r3 <- file.path(quantdir,"Nutl_t24_r3","quant.sf")

t48_r1 <- file.path(quantdir,"Nutl_t48_r1","quant.sf")
t48_r2 <- file.path(quantdir,"Nutl_t48_r2","quant.sf")
t48_r3 <- file.path(quantdir,"Nutl_t48_r3","quant.sf")


#merging the imported files into a single container (files) to be used to build the coldata
files <- c(t0_r1,t0_r2,t0_r3,t8_r1,t8_r2,t8_r3,t12_r1,t12_r2,t12_r3,t16_r1,t16_r2,t16_r3,t24_r1,t24_r2,t24_r3,t48_r1,t48_r2,t48_r3)
#creating the names variable
names <- c ("Nutl_t0_r1","Nutl_t0_r2","Nutl_t0_r3","Nutl_t8_r1","Nutl_t8_r2","Nutl_t8_r3","Nutl_t12_r1","Nutl_t12_r2","Nutl_t12_r3","Nutl_t16_r1","Nutl_t16_r2","Nutl_t16_r3","Nutl_t24_r1","Nutl_t24_r2","Nutl_t24_r3","Nutl_t48_r1","Nutl_t48_r2","Nutl_t48_r3")
#creating the timepoint variable
timepoint <- c("t0","t0","t0","t8","t8","t8","t12","t12","t12","t16","t16","t16","t24","t24","t24","t48","t48","t48")
#creating replicate variable
replicate <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")

#building the coldata dataframe
coldata <- data.frame(files = files, names= names, timepoint=timepoint, replicate = replicate, stringsAsFactors=FALSE)
coldata

#activate tximeta and build the SummarizeExperiment (se) object
library(tximeta)
se <- tximeta(coldata)

#make the se (summarize experiment) ready for DGE analysis on Deseq2
gse <- summarizeToGene(se)
library(DESeq2)

#I can explore the gse object elements using:
head(rowRanges(gse)) #shows the genomic ranges
head(colData(gse)) #shows the information about the samples
head(assay(gse)) #shows the matrix of counts, abundance and lenght, as determined by salmon

#create a dds object starting from the summarized experiment (gse). In this case, I want the 
#samples to be sorted depending on the condition, specified when I created coldata (see, above).
#The gse object has the timepoints stored as characters (you can test this by gse$timepoint) 
#so Deseq2 needs to convert it into factors. This operation is done automatically (it will be prompted)
#and the success can be tested using dds$timepoint
dds <- DESeqDataSet(gse, design = ~ timepoint)

#prefiltering: removing rows that contain no counts to speed up computation
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
#in this specific case, I get from ∼40k rows to 18k
#additional filtering can be maintaining only those rows where at least 3 samples have at least
#a count of 10. The command is:
#keep <- rowSums(counts(dds)>=10)>=3

#variance stabilizing transformation and rlog
#see more details in http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#reading-in-data-with-tximeta
#in general vsd is better for larger datasets(n>30), rlog is better for smaller ones and when there is a wide range of sequencing depth
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
#blind=FALSE means that the experimental condition is not taken into account to compute the differences

#plotting the effect of the transformation
library("dplyr")
library("ggplot2")
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

norm_plot <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
norm_plot
ggsave("norm_plot.png", plot = norm_plot, device = "png", path = resdir, dpi = 300)
#plotted are the untransformed samples (log2), the vst correction and the rlog correction. 
#clearly, the rlog improves the linerity of the samples compared to the untransformed log2 data
#vst correction does not improve a lot the linearity on lower FC values (in this context)


## USING rlog for the following analysis ##

#assessing sample distances
#euclidean distance 
sampleDists <- dist(t(assay(rld)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$timepoint, vsd$replicate, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap_plot <-pheatmap(sampleDistMatrix,
                         clustering_distance_rows = sampleDists,
                         clustering_distance_cols = sampleDists,
                         col = colors)
pheatmap_plot
#exporting plot
ggsave("pheatmap.png", plot = pheatmap_plot, device = "png", path = resdir, dpi = 300)


#PCA. For other scaling, see the tutorial (link above, line 81 of the script)
plotPCA(vsd, intgroup = c("timepoint","replicate"))
pcaData <- plotPCA(rld, intgroup = c("timepoint","replicate"),returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = timepoint, shape = replicate)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rlog data")
PCA_plot
ggsave("PCA_plot.png", plot = PCA_plot, device = "png", path = resdir, dpi = 300)


#Differential Expression
dds <- DESeq(dds)
res <- results(dds)
res

#deseq computes the log fold change of replicate 1 vs replicate 2, as replicate was the 
#last argument when building the dds object (line 53)
#In general, the results for a comparison of any two levels of a variable can 
#be extracted using the contrast argument to results. 
#The user should specify three values: the name of the variable, the name of the level 
#for the numerator, and the name of the level for the denominator. 
#Here we extract results for the log2 of the fold change of t0 vs t12 and t0 vs t24
rest0t12 <- results(dds, contrast = c("timepoint", "t12","t0"))
rest0t12
summary(rest0t12)
plotMA(rest0t12, ylim=c(-6,6))

rest0t24 <- results(dds, contrast = c("timepoint", "t24","t0"))
rest0t24
summary(rest0t24)
plotMA(rest0t24)

rest0t12_stringent <- results(dds, alpha = 0.05, lfcThreshold = 1, contrast = c("timepoint","t12","t0"))
rest0t12_stringent
summary(rest0t12_stringent)
plotMA(rest0t12_stringent)

rest0t24_stringent <- results(dds, alpha = 0.05, lfcThreshold = 1, contrast = c("timepoint","t24","t0"))
rest0t24_stringent
summary(rest0t24_stringent)

rest12t24 <- results(dds, contrast = c("timepoint","t24","t12"))
rest12t24
summary(rest12t24) 

resFilt <- rest0t12[which(rest0t12$pvalue < 0.05 & abs(rest0t12$log2FoldChange)>1),]
write.csv(resFilt, file="rest0t12_filtered.csv")

#the standard parameters of results can be changed by:
#res.05 <- results(dds, alpha = 0.05)
#table(res.05$padj < 0.05)
#resLFC1 <- results(dds, lfcThreshold=1)
#table(resLFC1$padj < 0.1)
#which adjust the pvalue and the logFC to 0,05 and 1 respectively

#correct for multiple testing. Suppose res is the output of the results function on the dds object
#resSig <- subset(res, padj < 0.1)
#head(resSig[ order(resSig$log2FoldChange), ])
#head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

#Plotting
topGene <- rownames(rest0t12)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("timepoint"))


