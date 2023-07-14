library(rmarkdown)

getwd()
files <- c(
  "salmon_DESeq2/ZM_r1to6_analysis.Rmd",
  "salmon_DESeq2/Etop_r1to6_analysis.Rmd",
  "salmon_DESeq2/Noc_r1to6_analysis.Rmd",
  "salmon_DESeq2/Nutl_r1to6_analysis.Rmd",
  "salmon_DESeq2/DHCB_r1to6_analysis.Rmd",
  
)

output_format <- "html_document"


for (input_file in files) {
  render(input_file, output_format)
}
