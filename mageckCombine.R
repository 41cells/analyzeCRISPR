#! /opt/az/local/R/R-3.2.0/installdir/bin/Rscript

#############################################################################################
# Use: To aggregate FDR q-values from multiple mageck outputs into a single document
# Input: x.gene_summary.txt files produced from mageck 
# Output: Single output with combined q-values for postive and negative selection
# Requires: 
#############################################################################################

# Code example
# Rscript ~/bin/mageckCombine/mageckCombine.R log/nolog

print("mageckCombine.R is running")
args = commandArgs(trailingOnly = TRUE)
options(warn=-1)
logCall = as.character(args[1])

#####################
# A. Import dependencies
#####################
suppressPackageStartupMessages(library(stringr))
library(stringr)

#####################
# 1. Identify 'x.gene_summary.txt' files in the directory
#####################
allFiles = list.files()
allFiles = allFiles[grep('gene_summary.txt', allFiles)]

#####################
# 2. Build empty matrix based on first file dimensions
#####################

openFile = read.table(allFiles[1], header=TRUE, sep="\t", check.names=FALSE, as.is=T, stringsAsFactors=FALSE)
outputPos = matrix(nrow = nrow(openFile), ncol = length(allFiles))
outputNeg = outputPos

#####################
# 2. Loop over all files and extract positive and negative selection FDR q-values
#####################

for (i in 1:length(allFiles)){
	openFile = read.table(allFiles[i], header=TRUE, sep="\t", check.names=FALSE, as.is=T, stringsAsFactors=FALSE)
	openFile = openFile[order(openFile$id),]
	outputNeg[,i] = openFile$fdr.neg
	outputPos[,i] = openFile$fdr.pos
	allFiles[i] = gsub('.gene_summary.txt', '', allFiles[i])
}

rownames(outputNeg) = openFile$id
colnames(outputNeg) = allFiles
rownames(outputPos) = openFile$id
colnames(outputPos) = allFiles

#####################
# 3. Log transform q-values if 'log' parameter set
#####################

if (logCall == 'log'){
	outputPos = -log10(outputPos)
	outputNeg = -log10(outputNeg)
} else {
	outputPos = outputPos
	outputNeg = outputNeg
}

#####################
# 4. Generate output files
#####################
outputPos = data.frame(Gene = rownames(outputPos), outputPos)
outputNeg = data.frame(Gene = rownames(outputNeg), outputNeg)


dir.create('mageckCombine_output')
setwd(paste(getwd(), '/mageckCombine_output/', sep=''))
write.table(data.frame(allFiles), 'combined_files.txt', col.name=T, row.names=F, quote=FALSE,sep="\t")
write.table(outputPos, 'mageckCombine_positive_selection.txt', col.name=T,row.names=F, quote=FALSE,sep="\t")
write.table(outputNeg, 'mageckCombine_negative_selection.txt', col.name=T,row.names=F, quote=FALSE,sep="\t")
print("mageckCombine.R is complete")
print("Files are stored in mageckCombine_output")
