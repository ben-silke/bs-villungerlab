class RFileWriter():
    treatment: str
    directory: str
    file_location: str

    short_times = ("c(0, 8, 12, 16, 24, 48)", [8,12,16,24,48])
    long_times = ("c(0, 16, 20, 24, 36, 48)", [16,20,24,36,48])

    all_replicates: bool

    time_dict = {
        'ZM': long_times,
        'Noc': long_times,
        'DHCB': long_times,
        'Nutl': short_times,
        'Etop': short_times
    }

    _r = "{r}"
    _r_false_include = "{r include=FALSE}"
    _r_setup = "{r setup, include=FALSE}"


    def __init__(self, treatment: str, directory: str, file_location: str, all_replicates: bool) -> None:
        self.treatment = treatment
        self.directory = directory
        self.all_replicates = all_replicates
        self.file_location = file_location

    def write_markdown_file(self, include_data_create):
        outline = self.markdown_outline(include_data_create)
        tabset_detail = "{.tabset}"
        content = f"""
{outline}
# {self.treatment} Analysis {tabset_detail}
## Exploratory Analysis
{self.write_intro_analysis()}
"""

        for time in self.time_dict[self.treatment][1]:
            variable_name = f'results_t{0}_t{time}'

            time_content = f"""

            
## Time -- {self.treatment} T0 to T{time} {tabset_detail}
{self.exploration_analyis(0, time)}

{self.differential_expression_analysis(0,time)}

{self.write_batch_limma_analysis(variable_name, time, self.treatment)}

"""
            content = content + time_content
        
        # file = f'{self.directory}/{self.treatment}_analysis.Rmd'
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"

        file = f'{self.treatment}_{replicate_file_name}_analysis.Rmd'

        with open(file, 'w') as f:
            f.write(content)

    def write_intro_analysis(self):
        tabset_detail = "{.tabset}"
        variable = f'results_{self.treatment}'
        content = f"""

### PCA
#### VST - Variance Stablised Transformation
```{self._r}
vsd_{self.treatment} <- vst(dds, blind=FALSE)
plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))
```


#### rLog

```{self._r}
rld_{self.treatment} <- rlog(dds, blind=FALSE)
plotPCA(rld_{self.treatment}, intgroup=c('batch', 'timepoint'))
```


### Heatmap

#### Euclidian Distance

{self.write_heatmap(self.treatment)}

#### Poission Distance
{self.write_poi_heatmap(self.treatment)}



        
We can use the package limma to account for the batch effect:
https://bioconductor.org/packages/release/bioc/html/limma.html

We can see that there is a clear batch effect.
```{self._r}
vsd_{self.treatment} <- vst(dds, blind=FALSE)
plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))

```

Using limma we can normalise the counts for this.
```{self._r}
mat_{self.treatment} <- assay(vsd_{self.treatment})
mm_{self.treatment} <- model.matrix(~timepoint, colData(vsd_{self.treatment}))
mat_{self.treatment} <- limma::removeBatchEffect(mat_{self.treatment}, batch=vsd_{self.treatment}$batch, design=mm_{self.treatment})
assay(vsd_{self.treatment}) <- mat_{self.treatment}
plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))

```
    

"""
        return content
    
    def exploration_analyis(self, start, end):
        tabset_detail = "{.tabset}"
        timepoint_name = f'timepoint_t{end}_vs_t{start}'
        variable_name = f'results_t{end}_t{start}'

        content = f"""
### Exploratory Analysis {tabset_detail}

#### Summary
##### General Summary
```{self._r}
{variable_name} <- results(dds, contrast = c("timepoint", "t{end}", "t{start}"))

{variable_name} = add_annotations_to_results({variable_name})

summary({variable_name})

```
##### Restricted Values
```{self._r}
{variable_name}_restricted <-
results(
    dds,
    alpha = 0.1,
    ## We are looking for genes which have at least doubled, or decreased by more than half
    lfcThreshold = 1,
    contrast = c("timepoint", "t{end}", "t{start}")
)
summary({variable_name}_restricted)

```
{self.write_ma_plot("%s_restricted" % variable_name)}

"""
        return content

    def differential_expression_analysis(self, start, end):
        tabset_detail = "{.tabset}"
        timepoint_name = f'timepoint_t{end}_vs_t{start}'
        variable_name = f'results_t{end}_t{start}'

        content = f"""
### Differential Expression Analysis {tabset_detail}
It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

#### APEGLM
```{self._r}

{variable_name}_shrunk_apeglm <- lfcShrink(dds, coef="{timepoint_name}", type="apeglm")
summary({variable_name}_shrunk_apeglm)

```
###### MA - Plot
{self.write_ma_plot("%s_shrunk_apeglm" % variable_name)}

###### V-Plot
{self.write_vplot("%s_shrunk_apeglm" % variable_name)}


#### ASHR
```{self._r}
{variable_name}_shrunk_ashr <- lfcShrink(dds, coef="{timepoint_name}", type="ashr")
summary({variable_name}_shrunk_ashr)

```
###### MA - Plot
{self.write_ma_plot("%s_shrunk_ashr" % variable_name)}

###### V-Plot
{self.write_vplot("%s_shrunk_ashr" % variable_name)}

#### P-Adjusted Analysis
{self.write_p_adjusted_analysis(variable_name)}


"""
        return content


    def time_analysis(self, start, end):
        timepoint_name = f'timepoint_t{end}_vs_t{start}'
        variable_name = f'results_t{end}_t{start}'
        content = f"""
### Unrestricted
#### Summary
```{self._r}
{variable_name} <- results(dds, contrast = c("timepoint", "t{end}", "t{start}"))

{variable_name} = add_annotations_to_results({variable_name})

summary({variable_name})

```

### Restricted
#### Summary
```{self._r}
{variable_name}_restricted <-
results(
    dds,
    alpha = 0.1,
    ## We are looking for genes which have at least doubled, or decreased by more than half
    lfcThreshold = 1,
    contrast = c("timepoint", "t{end}", "t{start}")
)
summary({variable_name}_restricted)

```
{self.write_ma_plot("%s_restricted" % variable_name)}
{self.write_vplot("%s_restricted" % variable_name)}

### With LFC-Shrink
It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
#### APE GLM
```{self._r}

{variable_name}_shrunk_apeglm <- lfcShrink(dds, coef="{timepoint_name}", type="apeglm")
summary({variable_name}_shrunk_apeglm)

```

{self.write_ma_plot("%s_shrunk_apeglm" % variable_name)}
{self.write_vplot("%s_shrunk_apeglm" % variable_name)}

#### ashr
```{self._r}
{variable_name}_shrunk_ashr <- lfcShrink(dds, coef="{timepoint_name}", type="ashr")
summary({variable_name}_shrunk_ashr)

```

{self.write_ma_plot("%s_shrunk_ashr" % variable_name)}
{self.write_vplot("%s_shrunk_ashr" % variable_name)}

{self.write_p_adjusted_analysis(variable_name)}
        """

        return content
    
    
    def write_ma_plot(self, variable):
        content = """
MA Plot
```{r}

plotMA(%s)
```

        """ % variable

        return content

    def write_vplot(self, variable):
        content = f"""
V- Plot

```{self._r}
{variable} <- add_annotations_to_results({variable})
selected_genes <- as.character({variable}$symbol)

EnhancedVolcano({variable},
                lab = selected_genes,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8,8),
                ylab = expression(paste('-Log'[10],' adj P')),
                title = 'Differential expression for {variable}',
                pCutoff = 0.05,
                # at least double, or less than half
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)
```
"""
        return content

    def write_batch_limma_analysis(self, variable, timepoint, treatment):

        content = f"""

### Batch Effects 

##### Differential expression analysis with limma
After we have adjusted counts for the batch effect, 
we can see if different genes are highlighted for differential expression

```{self._r}
design_{variable} <- model.matrix(~ timepoint, colData(vsd_{treatment}))
fit_{variable} <- lmFit(mat_{self.treatment}, design_{variable})
# Replace 'treatment' and 'control' with your specific time points.
# Demonstrates difference between t0 and t{timepoint}
contrast.matrix <- makeContrasts(timepointt{timepoint}, levels=design_{variable})
fit2_{variable} <- contrasts.fit(fit_{variable}, contrast.matrix)
fit2_{variable} <- eBayes(fit2_{variable})

results_{variable} <- topTable(fit2_{variable}, adjust.method="BH", number=Inf)
results_{variable} <- add_annotations_to_results(results_{variable})

results_clean__{variable} <- results_{variable}[!is.na(results_{variable}$logFC), ]

results_ordered_{variable} <- results_clean__{variable}[order(results_clean__{variable}$logFC), ]
```

```{self._r}
selected_genes_{variable} <- as.character(results_ordered_{variable}$symbol)

EnhancedVolcano(
    results_ordered_{variable},
    lab = selected_genes_{variable},
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Batch controlled ZM t0_t{timepoint} (limma)',
    ylab = expression(paste('-Log'[10],' adj P')),       # Y-axis label
    # with adj p cutoff of 0.05
    pCutoff = 0.05,
    # at least double, or less than half
    FCcutoff = 1.0,
    pointSize = 3.0,
    labSize = 3.0
)

```

```{self._r_false_include}
results_ordered_{variable} <- add_annotations_to_results(results_ordered_{variable})
results_ordered_{variable}_df <- as.data.frame(results_ordered_{variable})
write.csv(results_ordered_{variable}_df, file = "../../../../results/batch_corrected_{variable}_data.csv")

```

"""
        return content


    def write_p_adjusted_analysis(self, variable):
        content = """
```{r}
resSig <- subset(%s, padj < 0.1) 
summary(resSig)
```

###### Strongest Down Regulated Genes
```{r}
head(resSig[ order(resSig$log2FoldChange), ])
```

###### Strongest up regulation
```{r}
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
```
        """ % variable

        return content



    def markdown_outline(self, include_data_create=False):
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"

        content = """
---
title: "Salmon Analysis - %s "
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")

load("../../../data/%s_data.RData")
dds
results <- results(dds)
resultsNames(dds)

add_annotations_to_results <- function(res) {
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                    keys=ens.str,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                    keys=ens.str,
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
return (res)
}
```
        """ % (f'{self.treatment}_{replicate_file_name}', f'{self.treatment}_{replicate_file_name}')

        if include_data_create:
            replicate = "1:6" if self.all_replicates else "1:3"
            times = self.time_dict[self.treatment][0]
            replicate_file_name = "r1to6" if self.all_replicates else "r1to3"
            function = """
add_annotations_to_results <- function(res) {
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                    keys=ens.str,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                    keys=ens.str,
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
return (res)
}
"""         
            included_content = 'r include_source, include=FALSE, code = readLines("r/src/utils.R"), code = readLines("r/src/pca_utils.R")'
            content = f"""
---
title: "Salmon Analysis - {self.treatment}: {replicate_file_name}"
output: html_document
---

```{self._r_setup}
knitr::opts_chunk$set(echo = TRUE)
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("limma")
library("RColorBrewer")
library("PoiClaClu")
library("vsn")
library("genefilter")

```


```{included_content}

times = {times}
treatment <- "{self.treatment}"
{self.file_location}

dds <- create_dds('{self.treatment}', data_directory, times, "salmon_quant", {replicate})
# Create the data and then save it
# save(dds, file = glue('r/data/', "{self.treatment}_{replicate_file_name}_data.RData"))

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results/{self.treatment}_{replicate_file_name}_data.csv")

# load("../../../data/%s_data.RData")
dds
results <- results(dds)
resultsNames(dds)

{function}

```
            """

        return content

    def write_heatmap(self, variable):

        content = f"""

```{self._r}
sampleDists_{variable} <- dist(t(assay(vsd_{variable})))

sampleDistMatrix_{variable} <- as.matrix( sampleDists_{variable} )
rownames(sampleDistMatrix_{variable}) <- paste( vsd_{variable}$timepoint, vsd_{variable}$batch, sep = "_" )
colnames(sampleDistMatrix_{variable}) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_{variable},
         clustering_distance_rows = sampleDists_{variable},
         clustering_distance_cols = sampleDists_{variable},
         col = colors)
```
"""
        return content


    def write_poi_heatmap(self, variable):
        content = f"""
```{self._r}
poisd_{variable} <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix_{variable} <- as.matrix( poisd_{variable}$dd )
rownames(samplePoisDistMatrix_{variable}) <- paste( vsd_{variable}$timepoint, vsd_{variable}$batch, sep = "_" )
colnames(samplePoisDistMatrix_{variable}) <- NULL
pheatmap(samplePoisDistMatrix_{variable},
         clustering_distance_rows = poisd_{variable}$dd,
         clustering_distance_cols = poisd_{variable}$dd,
         col = colors)
```
"""
        return content

class SalmonRFileWriter(RFileWriter):

    def write_r_file(self):
        replicate = "1:6" if self.all_replicates else "1:3"
        times = self.time_dict[self.treatment][0]
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"

        content = f"""
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")

times = {times}
treatment <- "{self.treatment}"
{self.file_location}

dds <- create_dds('{self.treatment}', data_directory, times, "salmon_quant", {replicate})
# Create the data and then save it
save(dds, file = glue('r/data/', "{self.treatment}_{replicate_file_name}_data.RData"))

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results/{self.treatment}_{replicate_file_name}_data.csv")
"""

        file = f"{self.treatment}_{replicate_file_name}_create_data.R"
        with open(file, 'w') as f:
            f.write(content) 

    
    

class StarRFileWriter(RFileWriter):
    def write_r_data_file(self):
        pass
    