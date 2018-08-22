#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
######################################
## download external data from GEO
##
#######################################
closeAllConnections()
require(htmltab)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nDownload SRA/fastq sequence data from GEO\n")
  cat("\tusage: downloadGEO.R /GEO/accession/number /path/to/output/folder/ [SE/PE]")
  cat("\t/GEO/accession/number: The accession number of GEO, such as GSE57872;\n")
  cat("\t/path/to/output/folder/: The out put folder where the fastq files will be saved\n")
  cat("\tSE/PE: (Optional) SE: single-end; PE: paired-end, default is 'SE'\n")
  cat("\teg:download.R GSE57872 /home/z5ouyang/GSE/\n")
  cat("\n\n")
  q()
}
strGSE <- args[1]
strOut <- paste(args[2],"/",sep="")
singleEnd <- TRUE
if(length(args)>2 && args[3]=="PE") singleEnd <- FALSE
if(substr(strGSE,1,3)!="GSE"){
  cat("Please provide GSE accession!\n")
  q()
}

setwd(strOut)
strWWW <- paste("https://ftp.ncbi.nlm.nih.gov/geo/series/",substr(strGSE,1,nchar(strGSE)-3),"nnn/",strGSE,"/soft/",strGSE,"_family.soft.gz",sep="")
system(paste("wget",strWWW))
system(paste("gunzip -f ",strGSE,"_family.soft.gz",sep=""))

A <- scan(paste(strGSE,"_family.soft",sep=""),character(),quiet=T,sep="\n")
gsmID <- gsub("\\^SAMPLE = ","",A[grep("\\^SAMPLE = ",A)])
sName <- gsub(" ","_",gsub("\\!Sample_title = ","",A[grep("\\!Sample_title = ",A)]))
sLink <- gsub("\\!Sample_relation = SRA: ","",A[grep("\\!Sample_relation = SRA: ",A)])
if(length(gsmID)!=length(sName) || length(gsmID)!=length(sLink)) stop("The number of data link are not match with the sample number.")

for(i in 1:length(sLink)){
  strSRA <- htmltab(sLink[i],1,1)[1,1]
  cat("Download",gsmID[i],sName[i],strSRA,"\n")
  strF <- paste(sName[i],"_",gsmID[i],".fastq.gz",sep="")
  if(file.exists(strF)) next
  if(!file.exists(paste(strSRA,".sra",sep="")))
    system(paste("wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr(strSRA,1,6),"/",
                 strSRA,"/",strSRA,".sra -O",strSRA,".sra",sep=""))
  strF <- paste(sName[i],"_",gsmID[i],".sra",sep="")
  system(paste("mv ",strSRA,".sra ",strF,sep=""))
  if(singleEnd) system(paste("fastq-dump --gzip",strF))
  else system(paste("fastq-dump -I --split-files --gzip",strF))
  system(paste("rm",strF))
}
cat("Download",strGSE,"successfully!\n")
