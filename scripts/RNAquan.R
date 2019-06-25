#!/usr/bin/env Rscript
############################################
## RNAquan.R
## 
#############################################
## the full path of each sample tag directory list
## Genome: mm10
## Results folder
rm(list=ls())
closeAllConnections()

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("Extract raw counts, TPM and alignment stats from tag directories les:\n")
  cat("Usage: RNAquan.R path/to/the/tagDir/list/file genome path/to/the/results/folder/\n")
  q()
}
strInput <- args[1]
strGenome <- "mm10"
if(length(args)>1) strGenome <- args[2]
strPath <- dirname(strInput)
if(length(args)>2) strPath <- args[3]
strOpt <- ""
if(length(args)>3) strOpt <- paste(args[4:length(args)],collapse=" ")

A <- read.table(strInput,as.is=T,sep="\t")#genes/exons -tbp 1 
strCmd <- paste("analyzeRepeats.pl rna ",strGenome," -condenseGenes -count exons -noadj ",strOpt," -d ",paste(A[,1],collapse = " ")," > ",strPath,"/HOMER.rawCount.txt",sep="")
system(strCmd)
strCmd <- paste("analyzeRepeats.pl rna ",strGenome," -count exons -tpm ",strOpt," -d ",paste(A[,1],collapse = " ")," > ",strPath,"/HOMER.rawTPM.txt",sep="")
system(strCmd)
######################################################
## matching TPM with condenseGene from raw counts
rawC <- read.table(paste(strPath,"/HOMER.rawCount.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
rawT <- read.table(paste(strPath,"/HOMER.rawTPM.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
rawT <- rawT[rownames(rawC),]
write.table(rawT,file=paste(strPath,"/HOMER.rawTPM.txt",sep=""),sep="\t",col.names = NA,quote=F)

######################################################
## stats
iStart <- 8
X <- read.table(paste(strPath,"/HOMER.rawCount.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
colnames(X) <- sapply(strsplit(colnames(X)," "),"[[",1)
stat.names <- c("starTotal","starUnique","starUniqueRate","starMulti","homerUniPos","homerTotal","Clonality","homerAvgLength","homerTotalCount","rc1","rc5","rc10","rc20")
stat <- matrix(NA,nrow=ncol(X)-iStart+1,ncol=length(stat.names),dimnames=list(colnames(X)[iStart:ncol(X)],stat.names))

for(i in colnames(X)[iStart:ncol(X)]){
cat(i,"\t")
  strLog <- list.files(i,"log$",full.names=T)
  if(length(strLog)==1){
    sStat <- scan(strLog,character(),sep="\t",quiet=T)
    sStat <- as.numeric(sStat[c(grep("input reads",sStat)+1,grep("Uniquely mapped reads number",sStat)+1,grep("Number of reads mapped to multiple loci",sStat)+1)])
    stat[i,1:4] <- c(sStat[1:2],sStat[2]/sStat[1],sStat[3])
  }
cat("tagInfo\t")
  res <- read.table(paste(i,"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T)
  stat[i,"homerUniPos"] <- res[1,"Unique.Positions"]#grepl("genome=",res[,1])
  stat[i,"homerTotal"] <- res[1,"Total.Tags"]#grepl("genome=",res[,1])
  stat[i,"Clonality"] <- as.numeric(gsub("averageTagsPerPosition=","",res[grepl("averageTagsPerPosition",res[,1]),1]))
  stat[i,"homerAvgLength"] <- as.numeric(gsub("averageTagLength=","",res[grepl("averageTagLength",res[,1]),1]))
  stat[i,"homerTotalCount"] <- sum(X[,i])
  stat[i,"rc1"] <- sum(X[,i]>0)
  stat[i,"rc5"] <- sum(X[,i]>4)
  stat[i,"rc10"] <- sum(X[,i]>9)
  stat[i,"rc20"] <- sum(X[,i]>19)
  cat("\n")
}
conn <- file(paste(strPath,"/HOMER.stats.txt",sep=""),"w")
cat("\t",sep="",file=conn)
write.table(stat,file=conn,sep="\t",quote=F)
close(conn)


########################################3
## from raw counts to TPM
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
if(F){
  X <- read.table(paste(strPath,"/HOMER.rawCount.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
  colnames(X) <- sapply(strsplit(colnames(X)," "),"[[",1)
  iStart <- 8
  
  l <- rep(0,iStart-1)
  for(i in colnames(X)[iStart:ncol(X)]){
    res <- read.table(paste(i,"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T)
    l <- c(l,as.numeric(gsub("averageTagLength=","",res[grepl("averageTagLength",res[,1]),1])))
  }
  
  for(i in iStart:ncol(X)){
    rate <- X[,i]/(X[,"Length"]-l[i]+1)
    ### the length of many exons with reads are shorter than the average length of reads 
    ### thus those ones TPM were set to 0
    rate[rate<=0] <- 0
    X[,i] <- rate*1e6/sum(rate)
  }
  colnames(X) <- c(colnames(X)[1:(iStart-1)],paste(colnames(X)[iStart:ncol(X)],"TPM"))
  conn <- file(paste(strPath,"/HOMER.rawTPM.txt",sep=""),"w")
  cat("\t",sep="",file=conn)
  write.table(X,file=conn,sep="\t",quote=F)
  close(conn)
}
