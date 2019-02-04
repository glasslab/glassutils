#!/bioinformatics/anaconda3_2019Jan/bin/Rscript
#############################
## scList.R
##
#############################

strPath <- "/data/singleCell/"
cat("\nThe available single cell data are:\n")
for(i in list.dirs(strPath,F,F)){
  for(j in list.dirs(paste(strPath,"/",i,sep=""),F,F)){
    cat(i,": ",j,"\n",sep="")
  }
}
cat("\n")
