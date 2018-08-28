#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
#########################################
## ATACmotif.R
##
#########################################
## findMotifsGenome.pl XXX_induced.txt genome ./XXX/ -size 200 -bg ./XXX_bg.txt
##
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nPlot the significant p-values of top known motifs of different conditions as a heatmap\n")
  cat("\tusage: ATACmotif.R /path/to/a/folder/contains/all/homer/motif/folders topN [-f /path/to/the/result/pdf -i /deNovo/motif/index]\n")
  cat("\t/path/to/a/folder/contains/all/homer/motif/folders: A path to a folder which contains homer motif analyses result folders,\n")
  cat("\t\twhich all different conditions would be plotted in a heatmap.\n")
  cat("\ttopN: A number indicate the top motifs to be selected.\n")
  cat("\t/path/to/the/result/pdf: (optional use -f) The pdf file where the result should be plotted, default, 'allMotif.pdf' in the folder as first specified parameter\n")
  cat("\t/deNovo/motif/index: (optional use -i) indicate the index of de novo motif to be plotted, seperated by ';', default: 1;2;3;...;topN\n")
  
  cat("\n\teg: ATACmotif.R /home/z5ouyang/motifAnalyses/ 5\n")
  cat("\n\n")
  q()
}
## input ----
strInput <- args[1]
topN <- as.numeric(args[2])
strPDF <- paste(strInput,"/allMotif.pdf",sep="")
deNovo <- 1:topN
require(getopt)
spec <- matrix(c("strPDF","f",2,"character",
                 "deNovo","i",2,"character"),byrow=TRUE, ncol=4)
options = getopt(spec,opt=commandArgs(trailingOnly=TRUE)[-c(1:2)])
if(sum(names(options)=="strPDF")==1) strPDF <- options$strPDF
if(sum(names(options)=="deNovo")==1) deNovo <- as.numeric(unlist(strsplit(options$deNovo,";")))

pdf(strPDF,width=9)
## known motif -------
motifP <- motifR <- c()
totalSeq <- c(target=0,bg=0)
for(i in list.dirs(strInput,recursive = F)){
  strMotif <- paste(i,"/knownResults.txt",sep="")
  if(file.exists(strMotif)){
    one <- read.table(strMotif,sep="\t",header=T,as.is=T,check.names=F,comment.char="")
    one <- one[1:min(topN,sum(one[,3]<0.01)),c(1,4,7,9),drop=F]
    dimnames(one) <- list(one[,1],c("motif",paste(basename(i),c("logP","target","bg"),sep="_")))
    # logP
    motifP <- merge(motifP,-one[,2,drop=F],by="row.names",all=T)
    rownames(motifP) <- motifP[,1]
    motifP <- motifP[,-1,drop=F]
    # enrichment
    one <- one[,3:4]
    motifR <- merge(motifR,matrix(as.numeric(gsub("%","",as.matrix(one))),ncol=ncol(one),dimnames=dimnames(one)),
                    by="row.names",all=T,sort=F)
    rownames(motifR) <- motifR[,1]
    motifR <- motifR[,-1,drop=F]
  }
}
# plot logP heatmap
if(length(motifP)==0) stop("Cannot locate any motif analyses result")
motifP[is.na(motifP)] <- 0
rownames(motifP) <- substr(rownames(motifP),1,apply(cbind(nchar(rownames(motifP)),15),1,min))
require(pheatmap)
require(RColorBrewer)
require(colorspace)
if(ncol(motifP)<=8){
  COL <- brewer.pal(n = 8, name ="Dark2")[1:ncol(motifP)]
}else{
  COL <-rainbow_hcl(ncol(motifP))
}
#COL <- brewer.pal(n = 8, name ="Dark2")[1:ncol(motifP)]
names(COL) <- colnames(motifP)

sMotif <- sort(motifP[motifP!=0])
heatCOL <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(min(length(sMotif)-2,10))
br <- c(0,sMotif[floor(seq(1,length(sMotif)-1,length.out=length(heatCOL)-1))]-min(diff(sMotif))/10,max(sMotif)+min(diff(sMotif))/10)
pheatmap(motifP,cluster_cols=F,annotation_colors=list(grp=COL),
         color = heatCOL,
         breaks=br,
         annotation_col=data.frame(row.names=colnames(motifP),
                                   grp=colnames(motifP)))
# plot bulble plot for enrichment
#print(motifR)
require(ggplot2)
X <- data.frame()
for(i in gsub("_bg","",grep("_bg$",colnames(motifR),value=T))){
  X <- rbind(X,data.frame(grp=i,motif=substr(rownames(motifR),1,apply(cbind(nchar(rownames(motifR)),15),1,min)),
                          enrichment=motifR[,paste(i,"target",sep="_")]/motifR[,paste(i,"bg",sep="_")],
                          bg=motifR[,paste(i,"bg",sep="_")]))
}

rownames(motifP) <- substr(rownames(motifP),1,apply(cbind(nchar(rownames(motifP)),15),1,min))
print(ggplot(X,aes(x=grp,y=motif))+geom_point(aes(size=enrichment,colour=bg))+scale_size_continuous(range = c(1,10))+
        scale_color_gradient(low="#fee5d9", high="#a50f15")+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line =element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.background = element_blank()))

## homer motif -------
require(htmltab)
ori <- par(mar=c(2,7,max(c(0,(30-length(deNovo)*2))),1)+0.1,mgp=c(0.5,0,0),tcl=-0.1)
for(i in list.dirs(strInput,recursive = F)){
  strMotif <- paste(i,"/homerResults.html",sep="")
  if(!file.exists(strMotif)) next
  # obtain the motif names
  motifID <- c()
  for(j in deNovo){
    res <- scan(paste(i,"/homerResults/motif",j,".info.html",sep=""),character(),quiet=T)
    motifID <- c(motifID,
                 paste(gsub("<H4>","",sapply(strsplit(head(res[grep("<H4>",res)],3),"\\/"),head,1)),collapse="/"))
  }
  # obtain the motif info
  motifs <- htmltab(strMotif,1,1)[rev(deNovo),]
  logP <- -as.numeric(motifs[,"log P-pvalue"])
  maxP <- ceiling(max(logP)/10)*10
  maxR <- maxP/3
  # significance p-value bar plot
  y <- barplot(logP,horiz=T,xlab="-log(p-value)",las=1,main=basename(i),xlim=c(0,maxP+maxR),col="#f7fcb9")
  text(rep(maxP/2,length(y)),y,rev(motifID))
  # target/bg ratio line plot
  COL <- c(frame="black",tg="#006d2c",bg="#74c476")
  axis(3,c(maxP,maxP+maxR/4,maxP+maxR/2,maxP+maxR*3/4,maxP+maxR),c("0%","25%","50%","75%","100%"),col=COL["frame"],col.axis=COL["frame"])
  mtext("Percentage in sequence",3,1,at=maxP+maxR,col=COL["frame"],adj=1)
  lines(as.numeric(gsub("%","",motifs[,"% of Targets"]))*maxR/100+maxP,y,col=COL["tg"])
  points(as.numeric(gsub("%","",motifs[,"% of Targets"]))*maxR/100+maxP,y,col=COL["tg"],pch=15)
  lines(as.numeric(gsub("%","",motifs[,"% of Background"]))*maxR/100+maxP,y,col=COL["bg"],lty=2)
  points(as.numeric(gsub("%","",motifs[,"% of Background"]))*maxR/100+maxP,y,col=COL["bg"],pch=16)
  lines(rep(maxP,2),c(0,max(y)*2),col=COL["frame"],lty=2)
  legend("bottomright",names(COL)[-1],col=COL[-1],lty=1:2,pch=15:16,box.lty=0,text.col=COL[-1])
}
par(ori)
dev.off()
cat("\nPlot motifs are successfully plotted\n")

