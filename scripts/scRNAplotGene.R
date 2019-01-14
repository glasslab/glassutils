#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
#############################
## scRNAplotGene.R
##
#############################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<4){
  cat("\n\nO'Young:\nPlot violin and tSNE of given gene expression for single cell RNAseq:\n")
  info <- data.frame(row.names=c("usage:","<assay type>:","<sc Name>:","<gene1>[,gene2;gene3,...]:","</path/to/a/.pdf>:","e.g.:"),
                     info=c("scRNAplotGene.R <assay type> <sc Name> <gene1>[,gene2,gene3,...] </path/to/a/.pdf>",
                            "The type of the assay, e.g. RNA",
                            "The name of the single cell name (detail from scList.R)",
                            "A list of gene symbols separated by ','",
                            "A pdf file name including the path",
                            "scRNAplotGene.R RNA hg38_GBM_Tumor1 SEPP1,CD48 scGENE.pdf"))
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
genes <- unlist(strsplit(args[3],","))
strPDF <- args[4]
if(substr(strPDF,nchar(strPDF)-2,nchar(strPDF))!="pdf") strPDF <- paste(strPDF,".pdf",sep="")
## process
X <- readRDS(strF[1])
allGenes <- rownames(X@data)
if(sum(!genes%in%allGenes)>0)
  cat("Following genes cannot be found for",args[1],args[2],":\n",
      paste(genes[!genes%in%allGenes],collapse="\n"),"\n")

## cluster identity
strF <- list.files(paste(strPath,args[1],args[2],sep="/"),"id$",full.names=T)[1]
if(!is.na(strF)){
  gID <- read.table(strF,sep="\t",as.is=T,header=F,row.names=1)
  X@ident <- plyr::mapvalues(x = X@ident, from = rownames(gID), to = gID[,1])
}

## plot
pdf(strPDF)#,width=4,height=4
for(i in genes[genes%in%allGenes]){
  print(VlnPlot(X, i, point.size.use=0.5,x.lab.rot=T))#,use.scaled=T,use.imputed=T,use.raw=T
  FeaturePlot(X, i, cols.use = c("#73737320", "#4a148680"),pt.size=0.2, 
              reduction.use ="tsne")
}
dev.off()
