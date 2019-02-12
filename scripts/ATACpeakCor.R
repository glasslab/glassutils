#!/usr/bin/env Rscript
####################################################
## ATACpeakCor.R
##
######################################################
## The input file includes the list of tag directories in three columns seperated by "/t". 
## First column: names of the merged peaks;
## Second column: Tag directories of replicates should be in the same row and sperated by ";"
## Third column: name of each replicates should be in the same row and sperated by ";"
## 
## The file used by RNAview can be used here
## 

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nPlot ATAC correlation among replications\n")
  cat("\tusage: ATACpeakCor.R /path/to/the/file/listed/sample/tag/directories genome [/path/to/the/result/folder/]\n")
  cat("\t/path/to/the/file/listed/sample/tag/directories: (An example attached at the end) The file contains path of tag directories, which organized as follow:\n")
  cat("\t\tAt least three columns which seperated by 'tab', columns after third will be ignored, each row is a condition, such as a cell type; or a treatment\n")
  cat("\t\tFirst column: the name of the condition, such as 'AMhigh_LPS03h' or 'Microglia_Tumor'\n")
  cat("\t\tSecond column: list full path of tag directories which belong to the condition specified in first column.\n")
  cat("\t\t\tThose tag directories are separated by ';'\n")
  cat("\t\tThird column: list the short name of each sample listed in the second column. Same order seperated by ';'\n")
  cat("\tgenome: mm10 or hg38.\n")
  cat("\t/path/to/the/result/folder/: (optional) A path to a folder where the out put correlation pdf files will be saved. Default the folder where the tag directory file is.\n")
  cat("\n\tAn example of the file (e.g. allTags.txt) contains tag directories:\n")
  cat("IMo_LPS06h\t/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_6h_B_ES_17_03_13_trimmed;/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_6h_A_ES_17_03_13_trimmed\tB;A\n")
  cat("IMo_LPS20h\t/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_20h_A_ES_17_03_13_trimmed;/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_20h_B_ES_17_03_13_trimmed\tA;B\n")
  cat("\n\teg: ATACpeakCor.R allTags.txt mm10 /home/z5ouyang/tagCor/\n")
  cat("\n\n")
  q()
}
require(gplots)
require(MASS)

strInput <- args[1]
strGenome <- args[2]
strOutput <- dirname(strInput)
if(length(args)>2) strOutput <- args[3]
strTmp <- paste(strOutput,"/ATACpeakCor_tmp",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)

## useful functions ------
limFun <- function(x){
  step <- diff(range(x))/20
  return(range(x)-c(step,-step))
}
plotColorTable <- function(X,col.scale=range(X,na.rm=TRUE),col.range=colorpanel(201,rgb(0.1,0.1,0.8),"white",rgb(0.8,0.1,0.1)),
                           col.mid=0,xlab="",ylab="",bText=TRUE,textCex=0.4,rCex=1,cCex=1,rLine=-2.8,cLine=-3,rClass=NULL,cClass=NULL,...){
  col.n <- length(col.range)
  X.c <- ncol(X)
  X.r <- nrow(X)
  step <- diff(sort(X))
  step <- min(step[step!=0])
  if(col.mid <= col.scale[1] || col.mid >= col.scale[2]){ 
    col.step <- seq(col.scale[1]-step/2,col.scale[2]+step/2,length.out=col.n+1)
  }else{
    a <- seq(col.scale[1]-step/2,col.mid-step/2,length.out=ceiling(col.n/2))
    b <- seq(col.mid+step/2,col.scale[2]+step/2,length.out=floor(col.n/2)+1)
    col.step <- c(a,b)
  }
  X.col <- matrix("white",nrow=nrow(X),ncol=ncol(X))
  for(i in 1:col.n)
    X.col[which(X>col.step[i] & X<col.step[i+1],arr.ind=TRUE)] <- col.range[i]
  plot(c(),c(),type="n",xlim=c(0,X.c+1),ylim=c(0,X.r+1),axes=FALSE,xlab=xlab,ylab=ylab,...)
  axis(1,1:X.c,colnames(X),tick=FALSE,line=cLine,cex.axis=cCex,las=2)
  axis(2,1:X.r,rownames(X),tick=FALSE,line=rLine,las=1,cex.axis=rCex)
  for(i in 1:X.r){
    rect((1:X.c)-0.5,i-0.5,(1:X.c)+0.5,i+0.5,
         border="gray88",col=X.col[i,])
    if(bText) text(1:X.c,i,gsub("NA","",paste(X[i,])),cex=textCex)
  }
  if(!is.null(rClass)){
    for(i in (1:length(rClass))[diff(as.numeric(as.factor(rClass)))!=0]){
      lines(c(0.5,X.c+0.5),rep(i+0.5,2),lwd=2)
    }
  }
  if(!is.null(cClass)){
    for(i in (1:length(cClass))[diff(as.numeric(as.factor(cClass)))!=0]){
      lines(rep(i+0.5,2),c(0.5,X.r+0.5),lwd=2)
    }
  }
}

## processing ------
exps <- read.table(strInput,sep="\t",comment.char="",as.is=T)

for(i in 1:nrow(exps)){
  cat("\n\n\nSTART:",exps[i,1],"...............................\n\n\n\n",sep="")
  tags <- unlist(strsplit(exps[i,2],";"))
  if(length(tags)<2) next
  strF <- c()
  for(j in tags){
    if(nchar(j)<10) next
    strF <- c(strF,paste(strTmp,"/",basename(j),"_peaks.txt",sep=""))
    strCMD <- paste("findPeaks",j,"-style factor -minDist 200 -size 200 -o",tail(strF,1))
    cat(strCMD,"\n")
    system(strCMD)
    #if(!file.exists(tail(strF,1))) system(strCMD)
  }
  strMerge <- paste(strTmp,"/",exps[i,1],"_mergedPeaks.txt",sep="")
  strCMD <- paste("mergePeaks",paste(strF,collapse = " "),">",strMerge)
  system(strCMD)
  #if(!file.exists(strMerge)) system(strCMD)
  strAnno <- gsub("mergedPeaks.txt$","annotatePeaks.txt",strMerge)
  strCMD <- paste("annotatePeaks.pl",strMerge,strGenome,"-d",gsub(";"," ",exps[i,2]),">",strAnno)
  system(strCMD)
  #if(!file.exists(strAnno)) system(strCMD)
  X <- read.table(strAnno,sep="\t",header=T,as.is=T,row.names=1,comment.char = "",quote="",check.names=F)
  reads <- log2(X[,19:ncol(X)]+1)
  sID <- colnames(reads) <- unlist(strsplit(exps[i,3],";"))#basename(sapply(strsplit(colnames(reads)," "),head,1))

  strCor <- paste(strOutput,"/",exps[i,1],"_cor.pdf",sep="")
  pdf(strCor,width=4,height=4)
  par(cex=0.5,mgp=c(1.5,0.5,0),mar=c(2.5,2.5,0,0)+0.5)
  COR <- matrix(1,nrow=ncol(reads),ncol=ncol(reads),dimnames=list(sID,sID))
  for(j in 1:(ncol(reads)-1)){
    for(k in (j+1):ncol(reads)){
      index <- apply(reads[,c(j,k)],1,max)>0
      tmp <- reads[index,c(j,k)]
      #marks=paste(X[index,15],X[index,9],substr(gsub("[[:punct:] ]","",X[index,7]),1,3),sep="_")
      tryM <- try(f1 <- kde2d(tmp[,1],tmp[,2],n=500),silent=T)
      if(is.null(names(tryM))){
        plot(c(),c(),xlab=sID[j],ylab=sID[k],xlim=limFun(f1$x),ylim=limFun(f1$y))
        index <- rep(T,nrow(tmp))
      }else{
        imageCOL <- c("#FFFFFFFF",colorpanel(20,"#3300FFFF","#00FFFFFF","#CCFF00FF"),colorpanel(20,"#CCFF00FF","#FF9900FF","#AA0000FF"))
        image(f1,col=imageCOL,xlab=sID[j],ylab=sID[k],xlim=limFun(f1$x),ylim=limFun(f1$y))
        imageZero <- diff(range(f1$z))/length(imageCOL)
        index <- apply(tmp,1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
      }
      points(tmp[index,1],tmp[index,2],col=imageCOL[2],pch=20)
      lines(range(c(limFun(f1$x),limFun(f1$y))),range(c(limFun(f1$x),limFun(f1$y))),col="gray")
      COR[j,k] <- COR[k,j] <- cor(tmp[,1],tmp[,2])
      legend("topleft",paste("r=",round(COR[j,k],3)),box.lty=0)
      #plot(c(),c(),xlab=sID[j],ylab=sID[k],xlim=limFun(f1$x),ylim=limFun(f1$y))
      #text(tmp[index,1],tmp[index,2],marks[index],cex=0.2)
    }
  }
  plotColorTable(round(COR,3),textCex=0.8)
  dev.off()
}



