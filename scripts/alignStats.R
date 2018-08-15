#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
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
  cat("\tResults: will be send to the standard out, use '>' to capture the results in file\n")
  cat("\tUsage: alignStats.R /data/mm10/BMDM/ATAC > /home/z5ouyang/alignResults.txt\n")
  cat("\tUsage: alignStats.R /home/z5ouyang/tagList > /home/z5ouyang/alignResults.txt\n")
  cat("\n\n")
  q()
}
strInput <- args[1]

if(sum(dir(dirname(strInput))==basename(strInput))==0){
  strDir <- unlist(read.table(strInput,as.is=T))
}else if(sum(dir(dirname(strInput))==basename(strInput))==1){
  strDir <- list.dirs(strInput,recursive=F)[-1]
  #strDir <- strDir[strDir!=strInput]
}else{
  stop("unknown input!")
}
######################################################
## stats
stat.names <- c("alignerTotal","alignerUnique","alignerUniqueRate","alignerMulti","homerUniPos","homerTotal","tagPosition","homerAvgLength","mitoNum","mitoRate")
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
  stat[i,"homerAvgLength"] <- as.numeric(gsub("averageTagLength=","",res[grepl("averageTagLength",res[,1]),1]))
  index <- res[,1]=="chrM"
  if(sum(index)==1){
    stat[i,"mitoNum"] <- as.numeric(res[index,"Total.Tags"])
  }else if(sum(index)==0){
    stat[i,"mitoNum"] <- 0
  }else{
    stop("ERROR: more than one chrM")
  }
  stat[i,"mitoRate"] <- stat[i,"mitoNum"]/ stat[i,"homerTotal"]
}
cat("\t",paste(colnames(stat),collapse="\t"),"\n",sep="")
for(i in rownames(stat)) cat(i,"\t",paste(stat[i,],collapse="\t"),"\n",sep="")

#write.table(stat,file=paste(dirname(strInput),"/HOMER.stats.txt",sep=""),sep="\t",quote=F,col.names=NA)
