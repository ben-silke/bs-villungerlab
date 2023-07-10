
library(org.Hs.eg.db)

get_annotation <- function() {
    genenames <- mapIds(org.Hs.eg.db,keys = rownames(dds),column = "SYMBOL",keytype="ENSEMBL")
    annotation <- data.frame(gene_name=genenames,
                            row.names=rownames(dds),
                            stringsAsFactors = FALSE)

    # head(annotation)
    return(annotation)
}


annotation <- get_annotation()