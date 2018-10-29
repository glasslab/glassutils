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
COL <- NULL
if(!is.null(strInfo)){
  tagInfo <- read.table(strInfo,sep="\t",row.names=1,as.is = T,comment.char = "")
  pClass <- c()
  for(i in rownames(tagInfo)){
    sID <- c(sID,paste(i,unlist(strsplit(tagInfo[i,2],";")),sep="_"))
    pClass <- c(pClass,rep(i,length(sID)-length(pClass)))
    if(ncol(tagInfo)>2) COL <- c(COL,tagInfo[i,3])
  }
  if(ncol(tagInfo)>2){
    names(COL) <- unique(pClass)
    print(COL)
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
write.table(cbind(idInUse=sID,tagName=basename(sapply(strsplit(colnames(rawC)," "),head,1)),group=pClass),file=paste(strOutput,"sID.txt",sep=""),sep="\t",row.names=F)
colnames(rawC) <- sID

distalC <- rawC[(is.na(rawTags$'Distance to TSS')|rawTags$'Distance to TSS'>distal)&grepl("Intergenic",rawTags$Annotation,ignore.case=T),]
if(distal<1000) distalC <- rawC[is.na(rawTags$'Distance to TSS')|rawTags$'Distance to TSS'>distal,]
distalC <- distalC[apply(distalC,1,function(x){return(sum(x>8))})>1,]
## DESeq2 -----
require(DESeq2)
pheno <- data.frame(row.names = colnames(distalC),grp=pClass)
D <- DESeqDataSetFromMatrix(countData=matrix(as.integer(distalC),nrow=nrow(distalC),dimnames=dimnames(distalC)),
                            colData=pheno,
                            design=as.formula(paste("~",paste(colnames(pheno),collapse="+"))))
dds <- DESeq(D,betaPrior=TRUE)
## extract target peaks ------
peakDef <- rawTags[rownames(distalC),1:4]
normP <- log2(counts(dds,normalized=T)+1)
allDBP <- c()
require(gplots)
require(MASS)
require(RColorBrewer)
require(colorspace)
if(is.null(COL)){
  if(length(unique(pClass))<=8){
    COL <- brewer.pal(n = 8, name ="Dark2")[1:length(unique(pClass))]
  }else{
    COL <-rainbow_hcl(length(unique(pClass)))
  }
  names(COL) <- unique(pClass)
}
pdf(paste(strOutput,"/pairwised_DCA.pdf",sep=""),width=4,height=4)
par(mar=c(2,2,0,0)+0.2,mgp=c(1,0.1,0),tcl=-0.05)
imageCOL <- c("#FFFFFFFF",colorpanel(20,"#3300FFFF","#00FFFFFF","#CCFF00FF"),colorpanel(20,"#CCFF00FF","#FF9900FF","#AA0000FF"))
for(i in unique(pClass)){
  peakID <- c()
  for(j in unique(pClass)){
    if(i==j) next
    res <- results(dds,contrast = c("grp",i,j))
    write.table(cbind(data.frame(res),contrast=paste(i,j,sep="-")),file=paste(strOutput,"/DCA_",j,".vs.",i,".txt"),sep="\t",quote=F,col.names=NA)
    peakID <- c(peakID,rownames(res)[!is.na(res[,"padj"])&res[,"padj"]<0.05&res[,"log2FoldChange"]>logFC])
    ## pairwised ploting
    x <- apply(normP[,pClass==i,drop=F],1,mean)
    y <- apply(normP[,pClass==j,drop=F],1,mean)
    xlim <- range(x)
    ylim <- range(y)
    xlab <- paste(i,": log2 Mean normalized tag",sep="")
    ylab <- paste(j,": log2 Mean normalized tag",sep="")
    tryM <- try(f1 <- kde2d(x,y,n=200),silent=T)
    if(is.null(names(tryM))){
      plot(c(),c(),xlab=xlab,ylab=ylab)
      index <- rep(T,nrow(normP))
    }else{
      image(f1,col=imageCOL,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
      imageZero <- diff(range(f1$z))/length(imageCOL)
      index <- apply(cbind(x,y),1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
    }
    points(x[index],y[index],col=imageCOL[2],pch=20)
    lines(range(c(xlim,ylim)),range(c(xlim,ylim)),col="gray")
    lines(range(c(xlim,ylim)),range(c(xlim,ylim))+logFC,col="gray",lty=2)
    lines(range(c(xlim,ylim)),range(c(xlim,ylim))-logFC,col="gray",lty=2)
    legend("topleft",paste(j,": ",sum(!is.na(res[,"padj"])&res[,"padj"]<0.05&res[,"log2FoldChange"]< -logFC)," peaks",sep=""),
           box.lty=0,text.col=COL[j])
    legend("bottomright",paste(i,": ",sum(!is.na(res[,"padj"])&res[,"padj"]<0.05&res[,"log2FoldChange"]>logFC)," peaks",sep=""),
           box.lty=0,text.col=COL[i])
  }
  peakID <- unique(peakID)
  allDBP <- c(allDBP,peakID)
  if(length(peakID)>50){
    write.table(peakDef[peakID,],file=paste(strOutput,i,"_target.txt",sep=""),
                sep="\t",quote=F,col.names=F)
  }
  write.table(peakDef[!rownames(peakDef)%in%peakID,],
              file=paste(strOutput,i,"_bg.txt",sep=""),
              sep="\t",quote=F,col.names=F)
}
dev.off()
allDBP <- unique(allDBP)
write.table(peakDef[allDBP,],file=paste(strOutput,"/allDCA.txt",sep=""),sep="\t",quote=F,col.names=F)
require(pheatmap)
subDBP <- allDBP[sample(length(allDBP),min(10000,length(allDBP)))]
#save(normP,allDBP,subDBP,COL,pClass,file="test.RData")
pdf(paste(strOutput,"/allDCA.pdf",sep=""),height = 9)#,onefile=FALSE
pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),labels_row=rep("",length(subDBP)),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(20),
         annotation_col = data.frame(row.names=colnames(normP),
                                     grp=pClass))
pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),scale="row",labels_row=rep("",length(subDBP)),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(20),
         annotation_col = data.frame(row.names=colnames(normP),
                                     grp=pClass))
pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),scale="row",labels_row=rep("",length(subDBP)),
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(20),cluster_cols=F,
         annotation_col = data.frame(row.names=colnames(normP),
                                     grp=pClass))
dev.off()
cat("\n\nATACdiff is finished successfully!\n\n")


