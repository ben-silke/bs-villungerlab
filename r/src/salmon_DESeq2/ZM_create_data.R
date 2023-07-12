
library("vsn")
library("genefilter")
source("r/src/utils.R")
source("r/src/pca_utils.R")

treatment <- "ZM"

data_directory = file.path('/Volumes/bs_external/villunger', glue('organised/{treatment}/output_salmon'))

times = c(0, 16, 20, 24, 36, 48)

data_directory
dds <- create_dds('ZM', data_directory, times, "salmon_quant", 1:6)
# Create the data and then save it
save(dds, file = glue('r/data/', glue(treatment, "ZM_data.RData")))

dds
# Exporting results to csv
res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered <- add_annotations_to_results(resOrdered)

head(resOrdered)


resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results/ZM_data.csv")


# Alternate method to export results
# library("ReportingTools")
# htmlRep <- HTMLReport(shortName="report", title="My report",
#                       reportDirectory="./report")
# publish(resOrderedDF, htmlRep)
# url <- finish(htmlRep)
# browseURL(url)