#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
####################################################
## ATACdiff.R
##
######################################################
## 
## 
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nIdentify the distal differential chromatin accessibility (DCA) on provided raw tag counts in peaks\n")
  cat("\tusage: ATACdiff.R /path/to/the/annotated/raw/tag/file distal [-f /path/to/the/file/listed/sample/tag/directories -s /group/code/list -o /path/to/the/result/folder/ -l logFC]\n")
  cat("\t/path/to/the/annotated/raw/tag/file: A path to the annotated peaks with '-noadj -d' tag directories, the result 'allRawTags.txt' from 'ATACquan.R'\n")
  cat("\tdistal: A number indicate the distance in bp from TSS to be considering distal region. Only the peaks in the distal region will be considered for analyses\n")
  cat("\t/path/to/the/file/listed/sample/tag/directories: (optional use -f, one of -f/-s is required) The file contains path of tag directories,\n")
  cat("\t\t which used to generate raw tag annotatioin by 'ATACquan.R', where the group codes for all samples are indicated.\n")
  cat("\t/group/code/list: (optional use -s, one of -f/-s is required) <group code1>;<group code2>[;group code3;...], specify all samples included in the annotated raw tag file.\n")
  cat("\t/path/to/the/result/folder/: (optional use -o) A path to a folder where the out put files will be saved.\n\t\tDefault: the folder where the annotated raw tag file is.\n")
  cat("\tlogFC: (optional use -l) A number indicates the cut-off of logFC. Defaul 1 (2 Fold Change).\n")
  
  cat("\n\teg: ATACdiff.R allRawTags.txt 1000 -f quanTags.txt -o /home/z5ouyang/DCA/\n")
  cat("\n\n")
  q()
}
## input ----
strTags <- args[1]
distal <- as.numeric(args[2])
strOutput <- paste(dirname(strTags),"/",sep="")
logFC <- 1
sID <- strInfo <- pClass <- NULL
require(getopt)
spec <- matrix(c("strInfo","f",2,"character",
                 "pClass","s",2,"character",
                 "output","o",2,"character",
                 "logFC","l",2,"character"),byrow=TRUE, ncol=4)
options = getopt(spec,opt=commandArgs(trailingOnly=TRUE)[-c(1:2)])
if(sum(names(options)=="strInfo")==1) strInfo <- options$strInfo
if(sum(names(options)=="pClass")==1) pClass <- unlist(strsplit(options$pClass,";"))
if(sum(names(options)=="output")==1) strOutput <- options$output
if(sum(names(options)=="logFC")==1) logFC <- options$logFC
if(is.null(strInfo)&&is.null(pClass)){
  stop("ERROR: one of -f/-s has to be provided for DCA analysis\n")
}
## get the sample information ----
if(!is.null(strInfo)){
  tagInfo <- read.table(strInfo,sep="\t",row.names=1,as.is = T)
  pClass <- c()
  for(i in rownames(tagInfo)){
    sID <- c(sID,paste(i,unlist(strsplit(tagInfo[i,2],";")),sep="_"))
    pClass <- c(pClass,rep(i,length(sID)-length(pClass)))
  }
}
if(is.null(sID)){
  id <- 1
  sID <- paste(pClass,id,sep="_")
  while(sum(duplicated(sID))>0){
    id <- id +1
    sID[duplicated(sID)] <- paste(substr(sID[duplicated(sID)],1,nchar(sID[duplicated(sID)])-1),id,sep="")
  }
}
## process ------
rawTags <- read.table(strTags,as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="",check.names=F)
rawC <- as.matrix(rawTags[,-(1:18)])
if(ncol(rawC)!=length(pClass)) stop(paste("ERROR:",strTags,"(",ncol(rawC),") contains different sample number from the group codes(",length(pClass),") (eigher -f or -s)."))
write.table(cbind(idInUse=sID,tagName=basename(sapply(strsplit(colnames(rawC)," "),head,1)),group=pClass),file=paste(strOutput,"sID.txt"),sep="\t",row.names=F)
colnames(rawC) <- sID
distalC <- rawC[is.na(rawTags$'Distance to TSS')|rawTags$'Distance to TSS'>distal,]
distalC <- distalC[apply(distalC,1,function(x){return(sum(x>8))})>1,]
## DESeq2 -----
require(DESeq2)
pheno <- data.frame(row.names = colnames(distalC),grp=pClass)
D <- DESeqDataSetFromMatrix(countData=matrix(as.integer(distalC),nrow=nrow(distalC),dimnames=dimnames(distalC)),
                            colData=pheno,
                            design=as.formula(paste("~",paste(colnames(pheno),collapse="+"))))
dds <- DESeq(D,betaPrior=TRUE)
## extract induced peaks ------
peakDef <- rawTags[rownames(distalC),1:4]
normP <- log2(counts(dds,normalized=T)+1)
allDBP <- c()
for(i in unique(pClass)){
  peakID <- c()
  for(j in unique(pClass)){
    if(i==j) next
    res <- results(dds,contrast = c("grp",i,j))
    peakID <- c(peakID,rownames(res)[!is.na(res[,"padj"])&res[,"padj"]<0.05&res[,"log2FoldChange"]>logFC])
  }
  peakID <- unique(peakID)
  allDBP <- c(allDBP,peakID)
  if(length(peakID)>50){
    write.table(peakDef[peakID,],file=paste(strOutput,i,"_induced.txt",sep=""),
                sep="\t",quote=F,col.names=F)
  }
  write.table(peakDef[!rownames(peakDef)%in%peakID,],
              file=paste(strOutput,i,"_bg.txt",sep=""),
              sep="\t",quote=F,col.names=F)
}
allDBP <- unique(allDBP)
write.table(peakDef[allDBP,],file=paste(strOutput,"/allDCA.txt",sep=""),sep="\t",quote=F,col.names=F)
require(pheatmap)
require(RColorBrewer)
COL <- brewer.pal(n = 8, name ="Dark2")[1:length(unique(pClass))]
names(COL) <- unique(pClass)
pdf(paste(strOutput,"/allDCA.pdf",sep=""),height = 9)#,onefile=FALSE
pheatmap(normP[allDBP,],annotation_colors=list(grp=COL),labels_row=rep("",length(allDBP)),
         color = colorRampPalette(c("navy", "gray90", "firebrick4"))(50),
         annotation_col = data.frame(row.names=colnames(normP),
                                     grp=pClass))
pheatmap(normP[allDBP,],annotation_colors=list(grp=COL),scale="row",labels_row=rep("",length(allDBP)),
         color = colorRampPalette(c("navy", "gray90", "firebrick4"))(50),
         annotation_col = data.frame(row.names=colnames(normP),
                                     grp=pClass))
dev.off()
cat("\n\nATACdiff is finished successfully!\n\n")


