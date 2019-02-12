#!/usr/bin/env Rscript
#######################################
## ATACdegMotif.R
##
#######################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<5){
  cat("\n\nO'Young:\nFinding motif from peaks close to given genes.\n")
  cat("\tusage: ATACdegMotif.R /path/to/a/peak/file /path/to/a/gene/list/file genome dist /path/to/motif/result/folder/ [/path/to/the/background/peak/file]\n")
  cat("\t/path/to/a/peak/file: a file contains peaks, at least 5 columns peak file separated by 'tab', [peak name,chr,start,end,strand],\n")
  cat("\t/path/to/a/gene/list/file: a file contains gene definition (can be SYMBOL, NM_..., ENS..., accession number), one gene a row\n")
  cat("\tgenome: mm10 or hg38.\n")
  cat("\tdist: A number indicate the distance from given gene TSS where peaks are included for motif analyses.\n")
  cat("\t/path/to/motif/result/folder/: The folder where the result files should be saved\n")

  cat("\n\teg: ATACdegMotif.R IMo_LPS06h_induced.txt DEGlist.txt mm10 30000 /home/z5ouyang/DEGmotif/\n")
  cat("\n\n")
  q()
}
strPeaks <- args[1]
strGene <- args[2]
strGenome <- args[3]
distal <- as.numeric(args[4])
strOut <- args[5]
strBG <- NULL
if(length(args)>5) strBG <- args[6]

## loading gene location information ----
if(grepl("mm",strGenome)){
  gID <- read.table("/bioinformatics/bioinformatics/homer410/data/accession/mouse2gene.tsv",as.is=T,sep="\t",quote="")
}else if(grepl("hg",strGenome)){
  gID <- read.table("/bioinformatics/bioinformatics/homer410/data/accession/human2gene.tsv",as.is=T,sep="\t",quote="")
}else{
  stop("Unknown genome!")
}
gLoc <- read.table(paste("/bioinformatics/bioinformatics/homer410/data/genomes/",strGenome,"/",strGenome,".tss",sep=""),as.is=T,sep="\t",quote="",row.names=1)

## get gene location ------------------
genes <- unlist(read.table(strGene,as.is=T))
genesExt <- gID[gID[,4]%in%gID[toupper(gID[,1])%in%toupper(genes),4],1] ## match genes in all different form
gLoc <- gLoc[rownames(gLoc)%in%genesExt,]
gLoc[,2] <- apply(gLoc[,2:3],1,mean) - distal
gLoc[,3] <- apply(gLoc[,2:3],1,mean) + distal

## select distal peaks
peaks <- read.table(strPeaks,sep="\t",as.is=T,row.names=1)
if(!dir.exists(strOut)) dir.create(strOut)
strLoc <- paste(strOut,"/gLoc.bed",sep="")
strPeak <- paste(strOut,"/peaks.bed",sep="")
write.table(gLoc[,1:3],file=strLoc,col.names = F,quote=F,row.names=F,sep="\t")
write.table(cbind(peaks[,1:3],rownames(peaks),rep(0,nrow(peaks)),peaks[,4]),file=strPeak,col.names = F,quote=F,row.names=F,sep="\t")
strGenePeak <- paste(strOut,"/genePeaks.bed",sep="")
system(paste("bedtools intersect -wa -u -a",strPeak,"-b",strLoc,">",strGenePeak))
nPeaks <- as.numeric(sapply(strsplit(system(paste("wc -l",strGenePeak),T)," "),head,1))
if(nPeaks<30){
  stop(paste("Only",nPeaks,"peaks overlap. At least 30 is required for the motif calling."))
}
if(is.null(strBG)){
  system(paste("findMotifsGenome.pl",strGenePeak,strGenome,strOut,"-size given "))
}else{
  system(paste("findMotifsGenome.pl",strGenePeak,strGenome,strOut,"-size given -bg",strBG))
}
#-mknown /home/z5ouyang/src/data/all_threshold_0.5.motif





