# 
# 
# # Requirements to use RSubread - featureCounts to read in STAR data
# 
# 1. Aligned reads
# 2. genomic annotated features in GTF or GFF

library(Rsubread)


###################################################
### code chunk number 4: Rsubread.Rnw:125-135
###################################################
ann <- data.frame(
  GeneID=c("gene1","gene1","gene2","gene2"),
  Chr="chr_dummy",
  Start=c(100,1000,3000,5000),
  End=c(500,1800,4000,5500),
  Strand=c("+","+","-","-"),
  stringsAsFactors=FALSE)
ann
fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE