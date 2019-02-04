#!/bioinformatics/anaconda3_2019Jan/bin/Rscript
#############################
## scRNAmarker.R
##
#############################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3){
  cat("\n\nO'Young:\nPrint marker of a given cluster for single cell RNAseq:\n")
  info <- data.frame(row.names=c("usage:","<assay type>:","<sc Name>:","<cluster name>:","[bg cluster name]","e.g.:"),
                     info=c("scRNAmarker.R <assay type> <sc Name> <cluster name> [bg cluster name]",
                            "The type of the assay, e.g. RNA",
                            "The name of the single cell name (detail from scList.R)",
                            "A cluster identity name (detail from scRNAtsne.R)",
                            "(optional) A backgroup cluster identity name (detail from scRNAtsne.R), default all others",
                            "scRNAmarker.R RNA hg38_GBM_Tumor1 DAMs > DAMs_markers.txt"))
  colnames(info) <- NULL
  print(info,right=F)
  cat("\n\n")
  cat("If you think an update on identities are needed, please inform the authors.\n\n")
  q()
}
strPath <- "/data/singleCell/"
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"rds$",full.names=T)
if(length(strF)==0) stop("The data does NOT exisit!")
strID1 <- args[3]
strID2 <- NULL
if(length(args)>3) strID2 <- args[4]
require(Seurat)
require(dplyr)
## process
X <- readRDS(strF[1])
## cluster identity
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"id$",full.names=T)[1]
if(!is.na(strF)){
  gID <- read.table(strF,sep="\t",as.is=T,header=F,row.names=1)
  X@ident <- plyr::mapvalues(x = X@ident, from = rownames(gID), to = gID[,1])
}
## print markers
markers <- FindMarkers(object = X, ident.1 = strID1, ident.2=strID2, min.pct = 0.5)
print(markers)


