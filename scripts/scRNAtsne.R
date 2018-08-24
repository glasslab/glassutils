#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
#############################
## scRNAtsne.R
##
#############################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3){
  cat("\n\nO'Young:\nPlot tSNE of single cell RNAseq with cluster identification:\n")
  info <- data.frame(row.names=c("usage:","<assay type>:","<sc Name>:","</path/to/a/.pdf>:","e.g.:"),
                     info=c("scRNAtsne.R <assay type> <sc Name> </path/to/a/.pdf>",
                            "The type of the assay, e.g. RNA",
                            "The name of the single cell name (detail from scList.R)",
                            "A pdf file name including the path",
                            "scRNAtsne.R RNA hg38_GBM_Tumor1 scTSNEidentity.pdf"))
  colnames(info) <- NULL
  print(info,right=F)
  cat("\n\n")
  q()
}
strPath <- "/data/singleCell/"
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"rds$",full.names=T)
if(length(strF)==0) stop("The data does NOT exisit!")
require(Seurat)
require(dplyr)
strPDF <- args[3]
if(substr(strPDF,nchar(strPDF)-2,nchar(strPDF))!="pdf") strPDF <- paste(strPDF,".pdf",sep="")
## process
X <- readRDS(strF[1])

## cluster identity
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"id$",full.names=T)[1]
if(!is.na(strF)){
  gID <- read.table(strF,sep="\t",as.is=T,header=F,row.names=1)
  X@ident <- plyr::mapvalues(x = X@ident, from = rownames(gID), to = gID[,1])
}

## plot
pdf(strPDF,width=9)#
print(TSNEPlot(object = X,do.return = T, pt.size = 0.5,do.label =T,plot.title=paste(length(X@ident),"cells")))
print(TSNEPlot(object = X,do.return = T, pt.size = 0.5,do.label =F,plot.title=paste(length(X@ident),"cells")))
dev.off()

cat("\n\ntSNE is successfully generated!\n")