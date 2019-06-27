#!/usr/bin/env Rscript
####################################################
## ATACidr.R
##
## please add "/gpfs/data01/glasslab/home/jtao/software/anaconda3/bin"
## and "/home/jtao/code/seq_merge_pipe" into the head of $PATH
######################################################
## The input file includes the list of tag directories in two columns seperated by "/t". 
## First column: names of the merged peaks;
## Second column: Tag directories of replicates should be in the same row and sperated by ";"
## 
## The file used by RNAview can be used here
## 

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\nCalculate the IDR peaks for each condition between 2 replications\n")
  cat("\tusage: ATACidr.R /path/to/the/file/listed/sample/tag/directories [/path/to/the/result/folder/]\n")
  cat("\t/path/to/the/file/listed/sample/tag/directories: (An example attached at the end) The file contains path of tag directories, which organized as follow:\n")
  cat("\t\tTwo columns which seperated by 'tab', columns after second will be ignored, each row is a condition, such as a cell type; or a treatment\n")
  cat("\t\tFirst column: the name of the condition, such as 'AMhigh_LPS03h' or 'Microglia_Tumor'\n")
  cat("\t\tSecond column: list 2 paths of tag directories which belong to the condition specified in first column.\n")
  cat("\t\t\tThose 2 tag directories are separated by ';'\n")
  cat("\t/path/to/the/result/folder/: (optional) A path to a folder where the out put idr files will be saved.\n\t\tDefault: the folder where the tag directory file is.\n")
  cat("\n\tAn example of the file (e.g. idrTags.txt) contains tag directories:\n")
  cat("IMo_LPS06h\t/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_6h_B_ES_17_03_13_trimmed;/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_6h_A_ES_17_03_13_trimmed\n")
  cat("IMo_LPS20h\t/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_20h_A_ES_17_03_13_trimmed;/data/scratch/ouyang/Lung/LPS_IP_ATAC_adpTrim/tag_directories/C57BL6_IMo_ATAC_lysisTn5_LPS_20h_B_ES_17_03_13_trimmed\n")
  cat("\n\teg: ATACidr.R idrTags.txt /home/z5ouyang/idr/\n")
  cat("\n\n")
  q()
}
strInput <- args[1] #"/home/z5ouyang/Eniko_Lung/Lung3/tagDirTrimView"#
strOutput <- dirname(strInput)
if(length(args)>1) strOutput <- args[2]
strTmp <- paste(strOutput,"/ATACidr_tmp/",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)

## 
exps <- read.table(strInput,sep="\t",comment.char="",as.is=T)
#one.title <- one <- cutoff <- c()
peakN <- c()
#names(peakN) <- exps[,1]
for(i in 1:nrow(exps)){
  cat("\n\n\nSTART: ",exps[i,1],"...............................\n\n\n\n",sep="")
  strD <- unlist(strsplit(exps[i,2],";"))
  if(length(strD)!=2){
    cat("The number of tag directories is NOT 2! escaping",i,"\n")
    next
  }
  strF <- c()
  for(j in strD){
    strF <- c(strF,paste(strTmp,basename(j),"_peaks.txt",sep=""))
    strCMD <- paste("findPeaks",j,"-L 0 -C 0 -fdr 0.99 -minDist 200 -size 200 -o",tail(strF,1))
    system(strCMD)
  }
  strTmpF <- paste(strTmp,exps[i,1],"_idr.txt",sep="")
  strCMD <- paste("idr.R",paste(strF,collapse = " "),"-p",strTmpF)
  system(strCMD)
  strIDRpng <- paste(strTmpF,".out.png",sep="")
  file.copy(strIDRpng,paste(strOutput,"/",exps[i,1],"_idr.png",sep=""),overwrite=T)
  strIDR <- paste(strOutput,"/",exps[i,1],"_idr.txt",sep="")
  file.copy(strTmpF,strIDR,overwrite=T)
  peakN <- c(peakN,sapply(strsplit(system(paste("wc -l",strIDR),T)," "),head,1))
}
write.table(cbind(exps[,1],peakN),file=paste(strOutput,"/IDRpeakNumber.txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
cat("\n\nATACidr finished successfully!\n")


