# Boilerplate code which can be reused to make the analysis for salmon data

library(DESeq2)

analyze_timepoints <- function(timepoint1, timepoint2, dds) {
  # print headers
  cat(paste0("# ", timepoint1, " to ", timepoint2, "\n"))
  
  # Unrestricted
  ## Summary
  cat("## Unrestricted\n")
  cat("### Summary\n")
  
  results_unrestricted <-
    results(dds, contrast = c("timepoint", timepoint1, timepoint2))
  print(summary(results_unrestricted))
  
  ## MA Plot
  cat("### MA Plot\n")
  print(plotMA(results_unrestricted))
  
  # Restricted
  ## Summary
  cat("## Restricted\n")
  cat("### Summary\n")
  
  results_restricted <-
    results(
      dds,
      alpha = 0.05,
      lfcThreshold = 1,
      contrast = c("timepoint", timepoint1, timepoint2)
    )
  print(summary(results_restricted))
  
  ## MA Plot
  cat("### MA Plot\n")
  print(plotMA(results_restricted))
}