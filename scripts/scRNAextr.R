#!/usr/bin/env Rscript
## scRNAextr.R
##
#########################
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(Matrix)))) install.packages("Matrix",repos="https://cran.cnr.berkeley.edu/")
if(!require(Matrix)||!require(optparse))
  stop("R packages of Matrix and optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=getwd(), 
              help="Output directory path [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/scRNA/outs/...matrix/ [options]\n\tsuch as: %prog /home/z5ouyang/scratch/JOS_10X/LG5_Brodman10_Microglia/LG5/outs/filtered_feature_bc_matrix",
                          option_list=option_list,prog="scRNAextr.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/scRNA/outs/...matrix/ is required.\n\tsuch as: ", call.=FALSE)
}

strSC <- args[1]
opt = parse_args(opt_parser,args[-1])
####################################################
cID <- as.character(unlist(read.table(paste(strSC,"/barcodes.tsv.gz",sep=""))))
gene <- read.table(paste(strSC,"/features.tsv.gz",sep=""))
counts <- as.matrix(readMM(paste(strSC,"/matrix.mtx.gz",sep="")))

dimnames(counts) <- list(gene[,2],cID)
write.csv(counts,file=paste(opt$out,"/",gsub("_$","",gsub("\\/","_",strSC)),".csv",sep=""))

cat("scRNAextr.R finished successfully!\n")
