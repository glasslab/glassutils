#!/usr/bin/env Rscript
############################################
## alignStats.R
## 
#############################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\nExtract alignment stats from tag directories:\n")
  cat("\tUsage: alignStats.R path/to/the/tagDir[or path/to/the/tagDir/list/file]\n")
  cat("\tpath/to/the/tagDir: A path of folder which contains the tag directories\n")
  cat("\tpath/to/the/tagDir/list/file: A path of a flie listed path to tag directories, one tag directory in one line.\n")
  cat("\t\tCan be usefull when tag directories are not located in the same folder, such as different runs.\n")
  cat("\teg: alignStats.R /data/mm10/BMDM/ATAC > /home/z5ouyang/alignResults.txt\n")
  cat("\teg: alignStats.R /home/z5ouyang/tagList > /home/z5ouyang/alignResults.txt\n")
  cat("\tResults: will be send to the standard out, use '>' to capture the results in file\n")
  cat("\n\n")
  q()
}
strInput <- args[1]

if(file.exists(strInput) && !dir.exists(strInput)){
  strDir <- unlist(read.table(strInput,as.is=T))
}else if(dir.exists(strInput)){
  strDir <- list.dirs(strInput,recursive=F)#[-1]
  #strDir <- strDir[strDir!=strInput]
}else{
  stop("unknown input!")
}
#cat("Extract alignment stats from:\n",paste(strDir,collapse="\n"),"\n",sep="\t")
######################################################
## stats
stat.names <- c("alignerTotal","alignerUnique","alignerUniqueRate","alignerMulti","homerUniPos","homerTotal","tagPosition","FragLength","peakSize","homerAvgLength","mitoNum","mitoRate")
stat <- matrix(NA,nrow=length(strDir),ncol=length(stat.names),dimnames=list(strDir,stat.names))

for(i in strDir){
  #cat(i,"\n")
  strLog <- list.files(i,"star.log$",full.names=T)
  if(length(strLog)==1){
    sStat <- scan(strLog,character(),sep="\t",quiet=T)
    sStat <- as.numeric(sStat[c(grep("input reads",sStat)+1,grep("Uniquely mapped reads number",sStat)+1,grep("Number of reads mapped to multiple loci",sStat)+1)])
    stat[i,1:4] <- c(sStat[1:2],sStat[2]/sStat[1],sStat[3])
  }
  strLog <- list.files(i,"bowtie2.log$",full.names=T)
  if(length(strLog)==1){
    sStat <- scan(strLog,character(),sep=" ",quiet=T)
    sStat <- as.numeric(sStat[c(grep("reads;",sStat)-1,
                                grep("exactly",sStat)-3,
                                grep(">1",sStat)-3)])
    stat[i,1:4] <- c(sStat[1:2],sStat[2]/sStat[1],sStat[3])
  }
  res <- read.table(paste(i,"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T)
  stat[i,"homerUniPos"] <- as.numeric(res[1,"Unique.Positions"])#grepl("genome=",res[,1])
  stat[i,"homerTotal"] <- as.numeric(res[1,"Total.Tags"])#grepl("genome=",res[,1])
  stat[i,"tagPosition"] <- res[1,"Total.Tags"]/res[1,"Unique.Positions"]
  stat[i,"FragLength"] <- as.numeric(gsub("fragmentLengthEstimate=","",res[grepl("fragmentLengthEstimate",res[,1]),1]))
  stat[i,"peakSize"] <- as.numeric(gsub("peakSizeEstimate=","",res[grepl("peakSizeEstimate",res[,1]),1]))
  stat[i,"homerAvgLength"] <- as.numeric(gsub("averageTagLength=","",res[grepl("averageTagLength",res[,1]),1]))
  if(file.exists(paste(i,"/tagInfo_with_M.txt",sep=""))) res <- read.table(paste(i,"/tagInfo_with_M.txt",sep=""),sep="\t",as.is=T,header=T)
  index <- res[,1]=="chrM"
  if(sum(index)==1){
    stat[i,"mitoNum"] <- as.numeric(res[index,"Total.Tags"])
    stat[i,"mitoRate"] <- stat[i,"mitoNum"]/ as.numeric(res[1,"Total.Tags"])
  }else if(sum(index)==0){
    stat[i,"mitoRate"] <- stat[i,"mitoNum"] <- 0
  }else{
    stop("ERROR: more than one chrM")
  }
}
cat("\t",paste(colnames(stat),collapse="\t"),"\n",sep="")
for(i in rownames(stat)) cat(i,"\t",paste(stat[i,],collapse="\t"),"\n",sep="")

#write.table(stat,file=paste(dirname(strInput),"/HOMER.stats.txt",sep=""),sep="\t",quote=F,col.names=NA)
