### R code from vignette source 'Rsubread.Rnw'

###################################################
### code chunk number 1: Rsubread.Rnw:70-73
###################################################
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)


###################################################
### code chunk number 2: Rsubread.Rnw:90-92
###################################################
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align.stat <- align(index="reference_index",readfile1=reads,output_file="alignResults.BAM",phredOffset=64)


###################################################
### code chunk number 3: Rsubread.Rnw:99-103
###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")
align.stat2 <- align(index="reference_index",readfile1=reads1,readfile2=reads2,
output_file="alignResultsPE.BAM",phredOffset=64)


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


###################################################
### code chunk number 5: Rsubread.Rnw:142-144
###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",annot.ext=ann,isPairedEnd=TRUE)
fc_PE


###################################################
### code chunk number 6: Rsubread.Rnw:152-174
###################################################
if(grepl("linux", R.version$os)) {
  md5.zip <- "ffd5036b36e25e9b61efc412e71820dd"
  URL <- "https://shilab-bioinformatics.github.io/cellCounts-Example/cellCounts-Example.zip"
  temp.file <- tempfile()
  temp.dir <- tempdir()
  downloaded <- tryCatch({
      download.file(URL, destfile = temp.file)
      tools::md5sum(temp.file) %in% md5.zip
    },
    error = function(cond){
      return(FALSE)
    }
  )
  if(!downloaded) cat("Unable to download the file.\n") 
} else downloaded <- FALSE

if(downloaded){
  unzip(temp.file, exdir=paste0(temp.dir,"/cellCounts-Example"))
  library(Rsubread)
  buildindex(paste0(temp.dir,"/chr1"), 
    paste0(temp.dir,"/cellCounts-Example/hg38_chr1.fa.gz"))
}


###################################################
### code chunk number 7: Rsubread.Rnw:176-185
###################################################
if(downloaded){
  sample.sheet <- data.frame(
    BarcodeUMIFile = paste0(temp.dir,"/cellCounts-Example/reads_R1.fastq.gz"), 
    ReadFile = paste0(temp.dir,"/cellCounts-Example/reads_R2.fastq.gz"),
    SampleName="Example", stringsAsFactors=FALSE
  )
  counts <- cellCounts(paste0(temp.dir,"/chr1"), sample.sheet,  nthreads=1,
    input.mode="FASTQ", annot.inbuilt="hg38")
}


###################################################
### code chunk number 8: Rsubread.Rnw:187-188
###################################################
if(downloaded) print(counts$sample.info)


###################################################
### code chunk number 9: Rsubread.Rnw:190-191
###################################################
if(downloaded) print(dim(counts$counts$Example))


###################################################
### code chunk number 10: Rsubread.Rnw:211-213
###################################################
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]


###################################################
### code chunk number 11: Rsubread.Rnw:238-239
###################################################
propmapped("alignResults.BAM")

