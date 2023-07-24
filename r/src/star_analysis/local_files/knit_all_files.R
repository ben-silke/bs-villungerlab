library(rmarkdown)


files <- c(
#          "DHCB_r1to6_analysis.Rmd",
           "Etop_r1to6_analysis.Rmd",
           "Noc_r1to6_analysis.Rmd",
           "Nutl_r1to6_analysis.Rmd",
           "ZM_r1to6_analysis.Rmd" 
)

output_format <- "html_document"

for (input_file in files) {
  render(input_file, output_format)
}


