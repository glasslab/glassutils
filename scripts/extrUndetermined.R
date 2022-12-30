#!/usr/bin/env Rscript
#################################
##extrUndetermined.R
##
###################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\nextrUndetermined.R /path/to/undetermined....fastq.gz path/to/the/sample/index/file\n\n")
  cat("\tsample index file is the same as the Nextseq sample file which include two columns separated by ','\n")
  cat("\t\tFirst row is ignored, which in Nextseq sample sheets are \n\t\t\t[Data]\n")
  cat("\t\tFirst column contains the sample name; Second column contains the index barcode\n")
  q()
}
strUndetermined <- args[1]
strIndex <- args[2]
strPath <- dirname(strUndetermined)
strUnd <- gsub("\\.gz$","",strUndetermined)
if(!file.exists(strUnd)){
  cat("gunzip",strUndetermined,"...\n")
  system(paste("gunzip",strUndetermined))
}
A <- read.csv(strIndex,as.is=T,skip=1)
for(i in 1:nrow(A)){
  cat("working on",A[i,1],"\n")
  strF <- paste(A[i,1],".fastq",sep="")
  strGZ <- paste(strF,".gz",sep="")
  if(!file.exists(strGZ)){
    if(!file.exists(strF)) system(paste("echo \"..*", A[i,2],"\" | awk '{gsub(\"_\",\"\\\\_\",$0);$0=\"^@\"$0\".*?(\\\\n.*){3}\"}1' | pcregrep -oM -f - ",strUnd," > ",strF,sep=""))
    system(paste("gzip",strF))
  }
}


