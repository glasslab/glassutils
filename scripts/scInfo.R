#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
#############################
## scInfo.R
##
#############################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nDetail information of a single cell project:\n")
  info <- data.frame(row.names=c("usage:","<assay type>:","<sc Name>:","e.g.:"),
                     info=c("scInfo.R <assay type> <sc Name>",
                            "The type of the assay, e.g. RNA",
                            "The name of the single cell name (detail from scList.R)",
                            "scInfo.R RNA hg38_GBM_Tumor1"))

  colnames(info) <- NULL
  print(info,right=F)
  cat("\n\n")
  q()
}
strPath <- "/data/singleCell/"
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"info$",full.names=T)
cat("\n")
if(length(strF)>0) system(paste("cat",strF))
cat("\n")