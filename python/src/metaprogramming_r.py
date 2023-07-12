
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
        'DCHB': long_times,
        'Nutl': short_times,
        'Etop': short_times
    }

    _r = "{r}"

    def __init__(self, treatment: str, directory: str, file_location: str, all_replicates: bool) -> None:
        self.treatment = treatment
        self.directory = directory
        self.all_replicates = all_replicates
        self.file_location = file_location

    def write_r_file(self):
        replicate = "1:6" if self.all_replicates else "1:3"
        times = self.time_dict[self.treatment][0]


        content = f"""
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")
{self.file_location}

times = {times}
treatment <- "{self.treatment}"


dds <- create_dds('{self.treatment}', data_directory, times, "salmon_quant", {replicate})
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "{self.treatment}_data.RData")))

res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)


resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results/{self.treatment}_data.csv")

        """

        # file = f'{self.directory}/{self.treatment}_create_data.R'
        file = f"{self.treatment}_create_data.R"
        with open(file, 'w') as f:
            f.write(content) 

    def write_markdown_file(self):
        outline = self.markdown_outline()
        tabset_detail = "{.tabset}"
        content = f"""
{outline}
# {self.treatment} Analysis {tabset_detail}

"""

        for time in self.time_dict[self.treatment][1]:
            time_content = f"""

## Time -- {self.treatment} T0 to T{time}
{self.time_analysis(0,time)}

"""
            content = content + time_content
        
        # file = f'{self.directory}/{self.treatment}_analysis.Rmd'
        file = f'{self.treatment}_analysis.Rmd'

        with open(file, 'w') as f:
            f.write(content) 

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

### With LFC-Shrink
It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
#### APE GLM
```{self._r}

{variable_name}_shrunk <- lfcShrink(dds, coef="{timepoint_name}", type="apeglm")
summary({variable_name}_shrunk)

```

{self.write_ma_plot("%s_shrunk" % variable_name)}

#### ashr
```{self._r}
{variable_name}_shrunk <- lfcShrink(dds, coef="{timepoint_name}", type="ashr")
summary({variable_name}_shrunk)

```

{self.write_ma_plot("%s_shrunk" % variable_name)}

{self.write_p_adjusted_analysis(variable_name)}
        """

        return content
    
    
    def write_ma_plot(self, variable):
        content = """
#### MA Plot
```{r}

plotMA(%s)
```

        """ % variable

        return content

    def write_p_adjusted_analysis(self, variable):
        content = """
### P-Adjusted Analysis
```{r}
resSig <- subset(%s, padj < 0.1) 
summary(resSig)
```

### Strongest Down Regulated Genes
```{r}
head(resSig[ order(resSig$log2FoldChange), ])
```

### Strongest up regulation
```{r}
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
```
        """ % variable

        return content



    def markdown_outline(self):
        outline = """

---
title: "Salmon Analysis - %s "
output: html_document
---

```{r setup, include=FALSE}
    knitr::opts_chunk$set(echo = TRUE)
    library(DESeq2)
    library("apeglm")
    library("ashr")
    library(org.Hs.eg.db)

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
        """ % self.treatment
        return outline
