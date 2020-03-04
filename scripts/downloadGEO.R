#!/usr/bin/env Rscript
######################################
## download external data from GEO
##
#######################################
parallel <- F
closeAllConnections()
#require(htmltab)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nDownload SRA/fastq sequence data from GEO\n")
  cat("\tusage: downloadGEO.R /GEO/accession/number /path/to/output/folder/ [SE/PE]\n")
  cat("\t/GEO/accession/number: The accession number of GEO, such as GSE57872;\n")
  cat("\t/path/to/output/folder/: The out put folder where the fastq files will be saved\n")
  cat("\tSE/PE: (Optional) SE: single-end; PE: paired-end, default is 'SE'\n")
  cat("\teg:downloadGEO.R GSE57872 /home/z5ouyang/GSE/\n")
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
if(!dir.exists(strOut)&&!dir.create(strOut)) stop(paste("'",strOut,"' cannot be created",sep=""))
setwd(strOut)
strWWW <- paste("https://ftp.ncbi.nlm.nih.gov/geo/series/",substr(strGSE,1,nchar(strGSE)-3),"nnn/",strGSE,"/soft/",strGSE,"_family.soft.gz",sep="")
system(paste("wget",strWWW))
system(paste("gunzip -f ",strGSE,"_family.soft.gz",sep=""))

A <- scan(paste(strGSE,"_family.soft",sep=""),character(),quiet=T,sep="\n")
gsmID <- gsub("\\^SAMPLE = ","",A[grep("\\^SAMPLE = ",A)])
sName <- gsub(" ","_",gsub("\\!Sample_title = ","",A[grep("\\!Sample_title = ",A)]))
sLink <- gsub("\\!Sample_relation = SRA: ","",A[grep("\\!Sample_relation = SRA: ",A)])
if(length(gsmID)!=length(sName) || length(gsmID)!=length(sLink)) stop("The number of data link are not match with the sample number.")

unSucc <- c()
for(i in 1:length(sLink)){
  strSRA <- htmltab::htmltab(sLink[i],1,1)[1,1]
  unSucc <- c(unSucc,i)
  cat("\n\n\nDownloading",gsmID[i],sName[i],strSRA,"\n")
  strF <- paste(sName[i],"_",gsmID[i],".fastq.gz",sep="")
  if((singleEnd && file.exists(strF))||(!singleEnd && file.exists(paste(sName[i],"_",gsmID[i],"_1.fastq.gz",sep="")))){
    unSucc <- head(unSucc,-1)
    next
  }
  if(F){
    strSID <- paste(sName[i],"_",gsmID[i],".sra",sep="")
    if(!file.exists(strSID)){
      strSRAfile <- paste(strSRA,".sra",sep="")
      if(!file.exists(strSRAfile)){
        strLink1 <- strLink2 <- sraSize <- NULL
        strLink1 <- paste("ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr(strSRA,1,6),"/",strSRA,"/",strSRA,".sra",sep="")
        strLink2 <- htmltab::htmltab(paste("https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=",strSRA,sep=""),
                            "//table[@class='geo_zebra run-viewer-download']",1)[1,"Name"]
        if(RCurl::url.exists(strLink1)){
          strCMD <- paste("wget -c ",strLink1," -O ",strSRAfile,sep="")
        }else if(RCurl::url.exists(strLink2)){
          strCMD <- paste('curl -L "',strLink2,'" > ',strSRAfile,sep="")
          sraSize <- as.numeric(gsub(",","",gsub("Kb","",htmltab(paste("https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=",strSRA,sep=""),
                                                                 "//table[@class='geo_zebra run-viewer-download']",1)[1,"Size"])))*1000
        }else{
          stop("\t",strSRA,"cannot be found.")
        }
        cat(strCMD,"\n")
        system(strCMD)
        iter <- 0
        while(!is.null(sraSize)&&iter<10&&abs(file.size(strSRAfile)-sraSize)>1500){
          iter <- iter+1
          cat("\tTrying",iter,"times\n")
          system(strCMD)
        }
        if(!file.exists(strSRAfile)||
           file.size(strSRAfile)<1000||
           (!is.null(sraSize)&&abs(file.size(strSRAfile)-sraSize)>1500))
          stop(paste(strSRA,"file was NOT downloaded successfully!\n"))
      }
      system(paste("mv ",strSRA,".sra ",strSID,sep=""))
    }
    if(singleEnd) a <- try(system(paste("fastq-dump --gzip",strSID)))
    else a <- try(system(paste("fastq-dump -I --split-files --gzip",strSID)))
    if(a==0){
      system(paste("rm",strSID))
      unSucc <- head(unSucc,-1)
    }
  }else{
    strSID <- paste(sName[i],"_",gsmID[i],sep="")
    if(singleEnd) strCMD <- paste("fastq-dump --gzip",strSRA)
    else strCMD <- paste("fastq-dump -I --split-files --gzip",strSRA)
    if(parallel){
      strQsub <- paste(strSRA,".sh",sep="")
      cat("#!/bin/bash\n#$ -N ",strSRA,
          "\n#$ -wd ",getwd(),
          "\n#$ -pe node 1\n#$ -l h_rt=50:00:00\n#$ -o ",getwd(),"/",strSRA,
          ".log\n#$ -e ",getwd(),"/",strSRA,
          ".log\n#- End UGE embedded arguments\nsource /etc/profile.d/modules_bash.sh\nmodule load sratoolkit/2.8.1-3\n",
          strCMD,"\n",
          "for file in ",strSRA,"*; do mv $file ${file/",strSRA,"/",strSID,"}; done",
          sep="",file=strQsub)
      system(paste("qsub",strQsub))
      next
    }
    a <- iter <- 10
    while(a!=0&&iter>0){
      cat("\t",iter,"times left to try\n")
      a <- try(system(strCMD))
      iter <- iter -1
    }
    if(a==0) unSucc <- head(unSucc,-1)
    for(one in list.files(pattern=strSRA)) system(paste("mv",one,gsub(strSRA,strSID,one)))
  }
}
if(parallel){
  cat("All downloading submitted!\n")
  q()
}
if(length(unSucc)==0){
  cat("Download",strGSE,"successfully!\n")
}else{
  cat("The following were failed, please re-execute the command:\n")
  for(i in unSucc) cat(gsmID[i],sName[i],"\n")
}
