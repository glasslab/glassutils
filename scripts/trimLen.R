#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
##################
## trimLen.R
##
##################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\nTrim fastqs for a given length\n")
  cat("\tusage: trimLen.R /full/path/to/the/fastq/folder/ <length>\n")
  cat("\tlength: the number of leading bp to be kept\n")
  cat("\teg: trimLen.R /data/archive/18-07-09-glass/AN/18-07-09-glass_cleveland_atac/ 50")
  cat("\n\tResult: A folder named trim will be created in the folder and where all the trimed fastqs will be located\n")
  cat("\n\n")
  q()
}
strPath <- args[1]
strLength <- args[2]

strOut <- paste(strPath,"/trim/",sep="")
if(!dir.exists(strOut)) dir.create(strOut)
for(i in list.files(strPath,"fastq.gz$",full.names=T)){
  strF <- paste(strOut,gsub("fastq\\.gz","trim.fastq",basename(i)),sep="")
  system(paste("homerTools trim -len",strLength,i))
  system(paste("mv ",i,".trimmed ",strF,sep=""))
  system(paste("gzip",strF))
}
