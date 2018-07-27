#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
###############################################################
## ATACtrim.R
## O'Young
##
###########################################################
## provide the 

rm(list=ls())
graphics.off()
closeAllConnections()
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("Trim ATAC primer ATCTCCGAGCCCACGAGAC\n")
  cat("usage: /home/z5ouyang/src/ATACtrim.R /full/path/to/the/fastq/folder/\n")
  q()
}
strPath <- args[1]
if(length(args)>1){
  strRes <- args[2]
}else{
  strRes <- paste(strPath,"/trimmed/",sep="")
}
if(!dir.exists(strRes)) {dir.create(strRes)}

strCmd <- "homerTools trim -3 ATCTCCGAGCCCACGAGAC -mis 3 -minMatchLength 4 -min 20"#
L <- l.names <- c()
for(i in list.files(strPath,"gz$")){
  strF <- paste(strPath,"/",i,sep="")
  strR <- paste(strF,".trimmed",sep="")
  if(!file.exists(strR)){
    system(paste(strCmd,strF))
  }
  A <- read.table(gsub("trimmed","lengths",strR),sep="\t",header=T,comment.char = "",row.names=1,as.is=T)
  l.names <- c(l.names,sapply(strsplit(i,"_S"),"[[",1))
  L <- cbind(L,A[,"Fraction"])
  system(paste("gzip",strR))
  file.rename(paste(strR,"gz",sep="."),paste(strRes,gsub("_S","_trimmed_S",i),sep=""))
}
colnames(L) <- l.names
write.csv(L,file=paste(strRes,"trimmed.length.csv",sep=""),quote=F)

