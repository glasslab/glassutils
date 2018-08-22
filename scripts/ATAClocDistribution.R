#!/bioinformatics/bioinformatics/anaconda3/bin/Rscript
################################################
## ATAClocDistribution.R
##
#############################################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  cat("\n\nO'Young:\nCalculate ATAC peak distribution on normalized annotated peaks by the locations (tss-promoter, distal, no in peaks) \n")
  cat("\tusage: atacPeaksLoc.R /full/path/to/the/annotated/normalized/peak/file distal [/full/path/to/the/result/pdf/file]\n")
  cat("\t/full/path/to/the/annotated/normalized/peak/file: the file path to the annotated normalized peak file with '-d' all tag directories\n")
  cat("\tdistal: A number indicate the distance in bp from TSS to be considering distal region.\n")
  cat("\t/full/path/to/the/result/pdf/file: (optional) the file name where the results will be plot in. Defaul: 'ATAClocDist.pdf' in the same folder of annotated file\n")
  cat("eg: ATAClocDistribution.R ATAC_LPS00h_annotatePeaks.txt 1000\n")
  cat("\n\n")
  q()
}
strAnno <- args[1]
distal <- as.numeric(args[2])
strPDF <- paste(dirname(strAnno),"/ATAClocDist.pdf",sep="")
## processing ----
normTags <- read.table(strAnno,as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="",check.names=F)
normCounts <- normTags[,-(1:18)]
sID <- basename(sapply(strsplit(colnames(normCounts)," "),head,1))
names(sID) <- LETTERS[1:length(sID)]

STATs <- matrix(0,nrow=3,ncol=length(sID),dimnames=list(c("Promoter",paste("Distal ",distal/1000,"k",sep=""),"Not in Peaks"),names(sID)))
STATs[3,] <- 1-apply(normCounts,2,sum)/10000000
STATs[1,] <- apply(normCounts[!is.na(normTags$'Distance to TSS')&normTags$'Distance to TSS'<=distal,],2,sum)/10000000
STATs[2,] <- apply(normCounts[is.na(normTags$'Distance to TSS')|normTags$'Distance to TSS'>distal,],2,sum)/10000000
COL <- c("#ef8a62","#2166ac","#4d4d4d")
pdf(strPDF,width=9)
par(mar=c(12,2,0,0)+0.2,mgp=c(0.5,0,0),tcl=-0.03)
plot(c(),c(),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
legend("top",rownames(STATs),fill=COL,horiz=T)
legend("bottom",paste(names(sID),sID,sep=":"))
barplot(STATs,las=2,col=COL)
dev.off()
colnames(STATs) <- sID
write.csv(t(STATs),file=gsub("pdf$","csv",strPDF))
cat("Successfully plot the distribution.\n")



