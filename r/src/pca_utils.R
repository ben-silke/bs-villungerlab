
library(org.Hs.eg.db)

get_annotation <- function(dds) {
    row_names <- rownames(dds)
    genenames <- mapIds(org.Hs.eg.db,keys = row_names,column = "SYMBOL",keytype="ENSEMBL")
    annotation <- data.frame(gene_name=genenames,
                            row.names=row_names,
                            stringsAsFactors = FALSE)

    # head(annotation)
    return(annotation)
}

