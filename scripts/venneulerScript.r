#!/usr/bin/Rscript

# reads in a Homer venn file and produces a plot using the R venneuler package
# this script requires the R packages "venneuler" and "rJava"

### imports ###
library("rJava")
library("venneuler")
### end imports ###

### functions ###
# trims leading and trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

getVennWeights <- function(inputPath) {
# inputs: path to a Homer venn file
# outputs: a named vector accepted by the venneuler package venneuler method, text for a legend
	letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
	# read in put file
	data <- readLines(trim(inputPath))
	headerTokens <- strsplit(data[1], "\t")[[1]]
	header <- headerTokens[1:(length(headerTokens) - 2)]

	headerLength <- length(header)

	indexLetterDict <- letters[1:headerLength]
	
	legendText <- vector()
	for (i in 1:length(header)) {
		legendText <- c(legendText, paste(indexLetterDict[i],  " = ",  header[i], collapse=""))
	}

	weights <- vector()
	combinations <- vector()
	for (i in 2:length(data)) { # skip header
		tokens <- strsplit(data[i], "\t")[[1]]
		identityTokens <- tokens[1:headerLength]
		weights <- c(weights, as.numeric(tokens[headerLength + 1]))
		comb <- vector()
		for (j in 1:headerLength) {
			if (identityTokens[j] == "X")
				comb <- c(comb, indexLetterDict[j])
		}
		combinations <- c(combinations, paste(comb, collapse = "&"))
	}
	names(weights) <- combinations
	return(list("weights" = weights, "legendText" = legendText))
}

### end functions ###

# read in arguments
args <- commandArgs(trailingOnly = TRUE)

inputPath <- args[1]
outputPath <- args[2]

### HOW TO USE THE FUNCTION ###
results <- getVennWeights(inputPath) # call function
weights = results$weights # retrieve the weights
legendText= results$legendText # retrieve the legend text

print(legendText) # print the legend text to standard out

vd <- venneuler(weights) # compute venn diagram
jpeg(file=paste(trim(outputPath), "/vennDiagram.png", sep="", collapse="")) # sets file format and outputPath
plot(vd) # create plot at outputPath
pdf(file=paste(trim(outputPath), "/vennDiagram.pdf", sep="", collapse="")) # sets file format and outputPath
plot(vd) # create plot at outputPath
