#!/usr/bin/env Rscript
###################################
## idr.R
##
###################################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\nCalculate the IDR peaks for 2 replications\n")
  cat("\tusage: idr.R /path/to/peak1 /path/to/peak2 [-s idrThreshold] [-p idrThreshold]\n")
  cat("\te.g. idr.R rep1.peak rep2.peak\n")
  cat("\n\n")
  q()
}

threshold <- 0.05
strPeak1 <- args[1]
strPeak2 <- args[2]
strIDR <- ""
if(sum(grepl("-s",args))>0) threshold <- as.numeric(args[grep("-s",args)+1])
if(sum(grepl("-p",args))>0) strIDR <- args[grep("-p",args)+1]

peak1 <- read.table(strPeak1,sep="\t",as.is=T)
peak2 <- read.table(strPeak2,sep="\t",as.is=T)
strNarrow1 <- paste(strPeak1,"narrow",sep=".")
strNarrow2 <- paste(strPeak2,"narrow",sep=".")

## creat narrow peaks ----------------
write.table(cbind(peak1[,2:4],peak1[,1],rep(0,nrow(peak1)),peak1[,5],peak1[,8],
                  -log10(peak1[,12]),-log10(p.adjust(peak1[,12])),
                  #rep(-1,nrow(peak1)),rep(-1,nrow(peak1)),
                  round(apply(peak1[,3:4],1,diff)/2)),
            strNarrow1,sep="\t",col.names=F,row.names=F,quote=F)
write.table(cbind(peak2[,2:4],peak2[,1],rep(0,nrow(peak2)),peak2[,5],peak2[,8],
                  -log10(peak2[,12]),-log10(p.adjust(peak2[,12])),
                  #rep(-1,nrow(peak2)),rep(-1,nrow(peak2)),
                  round(apply(peak2[,3:4],1,diff)/2)),
            strNarrow2,sep="\t",col.names=F,row.names=F,quote=F)
## run idr ----------------------
if(nchar(strIDR)==0) strIDR <- paste(strPeak1,basename(strPeak2),"idr.peak",sep="_")
strOut <- paste(strIDR,"out",sep=".")
strCMD <- paste("idr --samples",strNarrow1,strNarrow2,
                "--output-file",strOut,"--plot",
                "--idr-threshold",threshold)
cat(strCMD,"\n")
system(strCMD)
## transfer idr peaks to homer peaks --------------------
idrPeak <- read.table(strOut,sep="\t",as.is=T)
idrPeak <- cbind(gsub(" ","",apply(idrPeak[,1:3],1,paste,collapse=":")),
                 idrPeak[,1:3],
                 idrPeak[,6:5],
                 idrPeak[,7])
colnames(idrPeak) <- c("#PeakID","chr","start","end","strand","idrScore","enrichScore")
conn <- file(strIDR,"w")
cat("#",strCMD,"\n",file=conn)
cat("# idr peak number:",nrow(idrPeak),"\n",file=conn)
cat("# Column Headers:\n",file=conn)
write.table(idrPeak,file=conn,row.names=F,quote=F,sep="\t")
close(conn)






















