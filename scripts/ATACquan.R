#!/bioinformatics/anaconda3_2019Jan/bin/Rscript
####################################################
## ATACquan.R
##
######################################################
## The input file includes the list of tag directories in two columns seperated by "/t". 
## First column: names of the merged peaks; where *_idrPeaks.txt can be located in the same location as the input file
## Second column: Tag directories of replicates should be in the same row and sperated by ";"
## 
## 
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nQuantify the tags in (provided) peaks (default the IDR peaks)\n")
  cat("\tusage: ATACquan.R /path/to/the/file/listed/sample/tag/directories genome [-f /path/to/peak/files -o /path/to/the/result/folder/ -d distal -s size]\n")
  cat("\t/path/to/the/file/listed/sample/tag/directories: The file contains path of tag directories.\n")
  cat("\tgenome: mm10 or hg38\n")
  cat("\t/path/to/peak/files: (optional use -f) path to the peak files (separated by ',', mergePeaks will be used) which will be used to quantify the ATAC tags.\n")
  cat("\t\tIf not provided, make sure all IDR peak files named as first column of tag directory list file plus '_idr.txt',\n")
  cat("\t\tand located in the same folder as tag directory list file.\n")
  cat("\t/path/to/the/result/folder/: (optional use -o) A path to a folder where the out put files will be saved.\n\t\tDefault: the folder where the tag directory list file is.\n")
  cat("\tdistal: (optional use -d) A number indicate the distance in bp from TSS to be considering distal region. defaul 1000\n")
  cat("\tsize: (optional use -s) A number indicate the Peak size. defaul given\n")

  cat("\n\teg: ATACquan.R quanTags.txt mm10 -o /home/z5ouyang/quan/\n")
  cat("\n\n")
  q()
}
## input ----
strTags <- args[1]
strGenome <- args[2]
tagInfo <- read.table(strTags,sep="\t",row.names=1,as.is = T)
strOutput <- strPath <- paste(dirname(strTags),"/",sep="")
strPeaks <- paste(strPath,rownames(tagInfo),"_idr.txt",sep="")
distal <- 1000
pSize <- "given"
require(getopt)
spec <- matrix(c("peaks","f",2,"character",
                 "output","o",2,"character",
                 "distal","d",2,"double",
                 "size","s",2,"double"),byrow=TRUE, ncol=4)
options = getopt(spec,opt=commandArgs(trailingOnly=TRUE)[-c(1:2)])
if(sum(names(options)=="peaks")==1) strPeaks <- unlist(strsplit(options$peaks,","))
if(sum(names(options)=="output")==1) strOutput <- options$output
if(sum(names(options)=="distal")==1) distal <- options$distal
if(sum(names(options)=="size")==1) pSize <- options$size

## get peak file (mergePeaks might be needed) ----
strTmp <- paste(strPath,"ATACquan_tmp/",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)
strPeak <- strPeaks
if(length(strPeaks)>1){
  strPeak <- paste(strTmp,"allPeaks",sep="")
  system(paste("mergePeaks -d",pSize,paste(strPeaks,collapse=" "),">",strPeak))#
}

## annotate peaks on tag directories ------
strCMD <- paste("annotatePeaks.pl",strPeak,strGenome,"-size",pSize,"-d",paste(unlist(strsplit(tagInfo[,1],";")),collapse=" "),">",paste(strPath,"allNormTags.txt",sep=""))
cat(strCMD,"\n")
system(strCMD)

strCMD <- paste("annotatePeaks.pl",strPeak,strGenome,"-size",pSize,"-noadj -d",paste(unlist(strsplit(tagInfo[,1],";")),collapse=" "),">",paste(strPath,"allRawTags.txt",sep=""))
cat(strCMD,"\n")
system(strCMD)

## plot the peak distrubution in location for each sample -----------
sID <- c()
for(i in rownames(tagInfo)) sID <- c(sID,paste(i,unlist(strsplit(tagInfo[i,2],";")),sep="_"))
normTags <- read.table(paste(strPath,"allNormTags.txt",sep=""),as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="")
normCounts <- normTags[,-(1:18)]
colnames(normCounts) <- sID
STATs <- matrix(0,nrow=3,ncol=length(sID),dimnames=list(c("Promoter",paste("Distal ",distal/1000,"k",sep=""),"Not in Peaks"),sID))
STATs[3,] <- 1-apply(normCounts,2,sum)/10000000
STATs[1,] <- apply(normCounts[!is.na(normTags$Distance.to.TSS)&normTags$Distance.to.TSS<=distal,],2,sum)/10000000
STATs[2,] <- apply(normCounts[is.na(normTags$Distance.to.TSS)|normTags$Distance.to.TSS>distal,],2,sum)/10000000
COL <- c("#ef8a62","#2166ac","#4d4d4d")
pdf(paste(strOutput,"ATACquan_dist.pdf",sep=""),width=9)
par(mar=c(12,2,0,0)+0.2,mgp=c(0.5,0,0),tcl=-0.03)
plot(c(),c(),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
legend("center",rownames(STATs),fill=COL)
barplot(STATs,las=2,col=COL)
dev.off()
cat("ATAC quantification is done successfully!\n")