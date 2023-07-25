class RFileWriter():
    treatment: str
    output: str
    file_location: str
    data_location: str

    short_times = ["c(0, 8, 12, 16, 24, 48)", [8,12,16,24,48], "c(8, 12, 16, 24, 48)"]
    long_times = ["c(0, 16, 20, 24, 36, 48)", [16,20,24,36,48], "c(16, 20, 24, 36, 48)"]

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


    def __init__(self, treatment: str, output: str = '', file_location: str = '', all_replicates: bool = False) -> None:
        self.treatment = treatment
        self.output = output
        self.all_replicates = all_replicates
        self.file_location = f"'{file_location}'"
        self.data_location = ''

    def write_markdown_file(self):
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"

        outline = self.markdown_outline()
        tabset_detail = "{.tabset}"
        self.replicate = replicate_file_name

        content = f"""
{outline}
# {self.treatment} Analysis {tabset_detail}
## Exploratory Analysis
{self.write_intro_analysis()}
"""

        for time in self.time_dict[self.treatment][1]:
            variable_name = f'results_{self.treatment}_t{0}_t{time}'

            time_content = f"""

            
## Time -- {self.treatment} T0 to T{time} {tabset_detail}
{self.exploration_analyis(0, time)}

{self.differential_expression_analysis(0,time)}

{self.write_batch_limma_analysis(variable_name, time, self.treatment)}

"""
            content = content + time_content
        

        file = f'{self.treatment}_{replicate_file_name}_analysis.Rmd'

        with open(file, 'w') as f:
            f.write(content)

    def write_intro_analysis(self):
        tabset_detail = "{.tabset}"
        variable = f'results_{self.treatment}_{self.replicate}'
        content = f"""

### PCA
#### VST - Variance Stablised Transformation
```{self._r}
vsd_{self.treatment} <- vst(dds_{self.treatment}_{self.replicate}, blind=FALSE)
plot <- plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))
plot + ggtitle("PCA plot highlighting batch differences for VST for {self.treatment}: {self.replicate}")
```


#### rLog

```{self._r}
rld_{self.treatment} <- rlog(dds_{self.treatment}_{self.replicate}, blind=FALSE)
plotPCA(rld_{self.treatment}, intgroup=c('batch', 'timepoint'))
```


### Heatmap

#### Euclidian Distance

{self.write_heatmap(self.treatment)}

#### Poission Distance
{self.write_poi_heatmap(self.treatment)}


### Batch Correction using _Limma_.
        
We can use the package _limma_ to account for the batch effect:
https://bioconductor.org/packages/release/bioc/html/limma.html

```{self._r}
vsd_{self.treatment} <- vst(dds_{self.treatment}_{self.replicate}, blind=FALSE)
plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))
```

```{self._r}
plot <- plotPCA(vsd_{self.treatment}, intgroup=c('batch'))
plot + ggtitle("PCA plot highlighting batch differences for VST for {self.treatment}: {self.replicate}")
```


Using limma we can normalise the counts for this.
```{self._r}
mat_{self.treatment} <- assay(vsd_{self.treatment})
mm_{self.treatment} <- model.matrix(~timepoint, colData(vsd_{self.treatment}))
mat_{self.treatment} <- limma::removeBatchEffect(mat_{self.treatment}, batch=vsd_{self.treatment}$batch, design=mm_{self.treatment})
assay(vsd_{self.treatment}) <- mat_{self.treatment}
plotPCA(vsd_{self.treatment}, intgroup=c('batch', 'timepoint'))

```

It is easier to view the change in timepoints here.
```{self._r}
plotPCA(vsd_{self.treatment}, intgroup=c('timepoint'))
```
    

"""
        return content
    
    def exploration_analyis(self, start, end):
        tabset_detail = "{.tabset}"
        timepoint_name = f'timepoint_t{end}_vs_t{start}'
        variable_name = f'results_{self.treatment}_t{end}_t{start}'

        content = f"""
### Exploratory Analysis {tabset_detail}

#### Summary
##### General Summary
```{self._r}
{variable_name} <- results(dds_{self.treatment}_{self.replicate}, contrast = c("timepoint", "t{end}", "t{start}"))

{variable_name} = add_annotations_to_results({variable_name})

summary({variable_name})

```
##### Restricted Values
```{self._r}
{variable_name}_restricted <-
results(
    dds_{self.treatment}_{self.replicate},
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
        variable_name = f'results_{self.treatment}_t{end}_t{start}'

        content = f"""
### Differential Expression Analysis {tabset_detail}
It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

#### APEGLM
```{self._r}

{variable_name}_shrunk_apeglm <- lfcShrink(dds_{self.treatment}_{self.replicate}, coef="{timepoint_name}", type="apeglm")
summary({variable_name}_shrunk_apeglm)

```
###### MA - Plot
{self.write_ma_plot("%s_shrunk_apeglm" % variable_name)}

###### V-Plot
{self.write_vplot("%s_shrunk_apeglm" % variable_name)}


#### ASHR
```{self._r}
{variable_name}_shrunk_ashr <- lfcShrink(dds_{self.treatment}_{self.replicate}, coef="{timepoint_name}", type="ashr")
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
        variable_name = f'results_{self.treatment}_t{end}_t{start}'
        content = f"""
### Unrestricted
#### Summary
```{self._r}
{variable_name} <- results(dds_{self.treatment}_{self.replicate}, contrast = c("timepoint", "t{end}", "t{start}"))

{variable_name} = add_annotations_to_results({variable_name})

summary({variable_name})

```

### Restricted
#### Summary
```{self._r}
{variable_name}_restricted <-
results(
    dds_{self.treatment}_{self.replicate},
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

{variable_name}_shrunk_apeglm <- lfcShrink(dds_{self.treatment}_{self.replicate}, coef="{timepoint_name}", type="apeglm")
summary({variable_name}_shrunk_apeglm)

```

{self.write_ma_plot("%s_shrunk_apeglm" % variable_name)}
{self.write_vplot("%s_shrunk_apeglm" % variable_name)}

#### ashr
```{self._r}
{variable_name}_shrunk_ashr <- lfcShrink(dds_{self.treatment}_{self.replicate}, coef="{timepoint_name}", type="ashr")
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

DESeq2::plotMA(%s, ylim=c(-3,3))
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
design_{self.treatment}_{timepoint} <- model.matrix(~ timepoint, colData(vsd_{treatment}))
fit_{self.treatment}_{timepoint} <- lmFit(mat_{self.treatment}, design_{self.treatment}_{timepoint})
# Replace 'treatment' and 'control' with your specific time points.
# Demonstrates difference between t0 and t{timepoint}
contrast.matrix <- makeContrasts(timepointt{timepoint}, levels=design_{self.treatment}_{timepoint})
fit2_{self.treatment}_{timepoint} <- contrasts.fit(fit_{self.treatment}_{timepoint}, contrast.matrix)
fit2_{self.treatment}_{timepoint} <- eBayes(fit2_{self.treatment}_{timepoint})

results_{self.treatment}_{timepoint} <- topTable(fit2_{self.treatment}_{timepoint}, adjust.method="BH", number=Inf)
results_{self.treatment}_{timepoint} <- add_annotations_to_results(results_{self.treatment}_{timepoint})

results_clean_{self.treatment}_{timepoint} <- results_{self.treatment}_{timepoint}[!is.na(results_{self.treatment}_{timepoint}$logFC), ]

results_ordered_{self.treatment}_{timepoint} <- results_clean_{self.treatment}_{timepoint}[order(results_clean_{self.treatment}_{timepoint}$logFC), ]
```

```{self._r}
selected_genes_{self.treatment}_{timepoint} <- as.character(results_ordered_{self.treatment}_{timepoint}$symbol)

EnhancedVolcano(
    results_ordered_{self.treatment}_{timepoint},
    lab = selected_genes_{self.treatment}_{timepoint},
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Batch controlled {self.treatment} t0_t{timepoint} (limma)',
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
# results_ordered_{self.treatment}_{timepoint} <- add_annotations_to_results(results_ordered_{self.treatment}_{timepoint})
# results_ordered_{self.treatment}_{timepoint}_df <- as.data.frame(results_ordered_{self.treatment}_{timepoint})
# write.csv(results_ordered_{self.treatment}_{timepoint}_df, file = "~/bs-villunger/results/batch_corrected_{self.treatment}_{timepoint}_data.csv")
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



    def markdown_outline(self):
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"
        variable = f"{self.treatment}_{replicate_file_name}"
        _r_setup = "{r setup, include=FALSE}"
        _functions = """
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
        content = f"""
---
title: "Salmon Analysis - {variable}"
output: html_document
---

```{_r_setup}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/bsilke/bs-villungerlab")
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library("glue")
library("Rsubread")
library("stringr")
library("DESeq2")
library(dplyr)
library(openxlsx)
load('~/bs-villungerlab/{self.output}/{self.treatment}_star_data.RData')

dds_{variable} <- ddseq_{self.treatment}
results <- results(dds_{variable})
resultsNames(dds_{variable})

{_functions}
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
poisd_{variable} <- PoissonDistance(t(counts(dds_{self.treatment}_{self.replicate})))

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
        self.replicate = replicate
        times = self.time_dict[self.treatment][0]
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"
        variable = f"{self.treatment}_{replicate_file_name}"

        self.data_location = f"{self.treatment}_{replicate_file_name}_salmon.RData"

        content = f"""
library("vsn")
library("genefilter")
# ("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

times = {times}
treatment <- "{self.treatment}"
data_directory = file.path({self.file_location})

{self.data_location}

dds_{variable} <- create_dds('{self.treatment}', data_directory, times, "salmon_quant", {replicate})
# Create the data and then save it
save(dds_{variable}, file = '{self.output}/{self.data_location}'))

res_{variable} <- results(dds_{variable})
resOrdered_{variable} <- res_{variable}[order(res_{variable}$padj),]
resOrdered_{variable} <- add_annotations_to_results(resOrdered_{variable})

head(resOrdered_{variable})

resOrderedDF_{variable} <- as.data.frame(resOrdered_{variable})
write.csv(resOrderedDF_{variable}, file = "results/{self.treatment}_{replicate_file_name}_data.csv")
"""

        file = f"{self.treatment}_{replicate_file_name}_create_data.R"
        with open(file, 'w') as f:
            f.write(content) 


class StarRFileWriter(RFileWriter):

    def write_r_file(self):
        replicate = "1:6" if self.all_replicates else "1:3"
        self.replicate = replicate
        times = self.time_dict[self.treatment][0]
        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"
        if not self.data_location:
            self.data_location = f"{self.treatment}_{replicate_file_name}_star.RData"
        for_loop = 'for (time in times) {'
        end_for_loop = '}'
        _time = '{time}'
        
        content = f"""
setwd("/Users/bsilke/bs-villungerlab")
library("glue")
library("stringr")
library("DESeq2")
library(DESeq2)
library("apeglm")
library("ashr")
library(EnhancedVolcano)
library(org.Hs.eg.db)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("limma")
library(dplyr)
library(openxlsx)

source("r/src/star_analysis/star_utils.R")
source("r/src/pca_utils.R")
source("r/src/utils.R")

times = {times}
treatment <- "{self.treatment}"
data_directory = file.path({self.file_location})
ddseq_{self.treatment} <- load_all_htseq_data(file.path(data_directory, 'all_{self.treatment}_htseq_encode_counts.tsv'))

# <- create_htseq_ddseq({self.treatment}, data_directory, times, {replicate})

save(ddseq_{self.treatment}, file = '{self.output}/{self.treatment}_star_data.Rdata')

{self.treatment}_workbook <- createWorkbook()
times = {self.time_dict[self.treatment][2]}
{for_loop}
    timepoint <- glue("timepoint_t{_time}_vs_t0")
    results_{self.treatment} <- lfcShrink(ddseq_{self.treatment}, coef=timepoint, type="apeglm")

    results_{self.treatment} <- subset(results_{self.treatment}, padj < 0.1)  # Restrict to values which are significant
    results_{self.treatment} <- add_annotations_to_results(results_{self.treatment})
    results_{self.treatment}_df <- as.data.frame(results_{self.treatment})
    addWorksheet({self.treatment}_workbook, glue("{self.treatment}_{_time}"))
    writeData({self.treatment}_workbook, glue("{self.treatment}_{_time}"), results_{self.treatment}_df, rowNames=TRUE)
{end_for_loop}

saveWorkbook({self.treatment}_workbook, "{self.output}/{self.treatment}_workbook.xlsx", overwrite = TRUE)

"""

        file = f"{self.treatment}_{replicate_file_name}_create_stardata.R"
        with open(file, 'w') as f:
            f.write(content) 



class ResultsSheetWriter:

    output: str
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


    def __init__(self, output: str, file_location: str, all_replicates: bool) -> None:
        self.output = output
        self.all_replicates = all_replicates
        self.file_location = file_location

        replicate_file_name = "r1to6" if self.all_replicates else "r1to3"
        self.replicate = replicate_file_name

    def write_result_creation_sheet(self):
        """
        Function to write the collection of data, for either star or for salmon data to write to spreadsheets
        
        """

        content = f"""

{self.write_intro()}
{self.write_create_data()}

#### TIMEPOINT t0-t16
{self.write_timepoint_workbook(16)}

#### TIMEPOINT t0-t24
{self.write_timepoint_workbook(24)}

#### TIMEPOINT t0-t48
{self.write_timepoint_workbook(48)}
"""

        file = f'all_results_{self.replicate}.R'

        with open(file, 'w') as f:
            f.write(content)

    def write_create_data(self):
        self.created_data = True
        raise NotImplementedError

    def write_timepoint_workbook(self, timepoint):
        timepoint_variable = f"timepoint_t{timepoint}_vs_t0"
        workbook_name = f"t0_t{timepoint}_workbook"
        treatments = [treatment for treatment in self.time_dict.keys() if timepoint in self.time_dict[treatment][1]]
        
        def write_treatment_timepoint_data(treatment, timepoint):
            

            content = f"""
addWorksheet({workbook_name}, "t0_t16_workbook")

### loop for each treatment Template:
results_{treatment}_{timepoint} <- lfcShrink(dds_{treatment}_{self.replicate}, coef="{timepoint_variable}", type="apeglm")
results_{treatment}_{timepoint} <- subset(results_{treatment}_{timepoint}, padj < 0.1)  # Restrict to values which are significant
results_{treatment}_{timepoint} <- add_annotations_to_results(results_{treatment}_{timepoint})
results_{treatment}_{timepoint}_df <- as.data.frame(results_{treatment}_{timepoint})

addWorksheet({workbook_name}, "results_apelgm_{treatment}_{timepoint}")
writeData({workbook_name}, "results_subset_dataframe", results_{treatment}_{timepoint}_df, rowNames=TRUE)
"""

            return content

        content = f"""
        
## TIMEPOINT t0-t{timepoint}- loop for each time stamp
{workbook_name} <- createWorkbook()

# Create results 
"""
        for treatment in treatments:
            content = content + write_treatment_timepoint_data(treatment, timepoint)


        content = content + f"""
saveWorkbook({workbook_name}, "{workbook_name}.xlsx", overwrite = TRUE)

"""
        return content

    def write_intro(self):

        content = f"""
library("vsn")
library("genefilter")
library(dplyr)
library(openxlsx)

setwd("/Users/bsilke/bs-villungerlab")
source("r/src/utils.R")
source("r/src/pca_utils.R")

short_times <- c(8,12,16,24,48)
long_times <- c(16,20,24,36,48)

# Create a named list
times <- list(
  'ZM' = long_times,
  'Noc' = long_times,
  'DHCB' = long_times,
  'Nutl' = short_times,
  'Etop' = short_times
)

treatments <- c("ZM", "DHCB", "Etop", "Noc", "Nutl")
path <- "{self.file_location}"
data_directory = file.path({self.file_location})

"""
        return content
    

class SalmonResultSheetWriter(ResultsSheetWriter):

    def write_create_data(self):

        def write_treatment_creation(treatment):
            content = f""" 

load("{self.output}/{treatment}_{self.replicate}_data.RData")
            
treatment = '{treatment}'
data_directory = file.path('{self.file_location}', glue('organised/{treatment}/output_salmon'))
dds_{treatment}_{self.replicate} <- create_dds(treatment, data_directory, times[treatment], "salmon_quant", 1:6)
save(dds_{treatment}_{self.replicate}, file = glue('{self.output}/', "{treatment}_{self.replicate}_data.RData"))

"""
            return content
        content = ""

        for treatment in self.time_dict.keys():
            content = content + write_treatment_creation(treatment)

        return content
