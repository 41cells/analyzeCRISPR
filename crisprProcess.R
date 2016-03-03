#! /opt/az/local/R/R-3.2.0/installdir/bin/Rscript

# Objective: To take counts file from CRISPR screening (Sabatini lab) and process for downstream workflows
# Input: Counts matrix from CRISPR screening with a single early time point
# Output: Normalized CRISPR ratios on the sgRNA level, gene level, mageck inputs

# Code structure: Rscript crisprProcess.R COUNTS_FILE INPUT_COLUMN_NAME NORMALIZATION(quantile/zmad) OUTPUT_DIRECTORY
# Code example: Rscript /gpfs/users/krxw569/bin/crisprProcess/crisprProcess.R breast.output.counts.txt

######################
# 1. Start-up
#####################
print("crisprProcess.R is running")
print("Checking options")

# Command line arguments into R
options(warn=-1)
args = commandArgs(trailingOnly = TRUE)

# Check if input arguments are correct
if (length(args) < 4) {
	stop('Requires 3 trailing options: Counts file, normalization method, and output directory')
} else {
	args = args
	print("Options are OK")
}

# Store object names
counts = as.character(args[1]); inputColumn = toupper(as.character(args[2])); method = as.character(args[3]); outputDir = as.character(args[4])

# Generate output directory
if (file.exists(outputDir)){
	outputDir = outputDir
} else {
	dir.create(outputDir)
}
print(paste('Output to be stored in ', outputDir, sep=""))


######################
# 2. Load global dependencies
#####################

# Libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(reshape2))
library(gplots)
library(stringr)
library(reshape2)

######################
# 3. Load counts file 
#####################

countsFile = read.table(counts, header=TRUE, as.is=TRUE, row.names=1, sep="\t")
setwd(outputDir)
colnames(countsFile) = toupper(colnames(countsFile))
colnames(countsFile) = gsub('\\.', '', colnames(countsFile))

######################
# 4. Extract gene names from dataset and place in first column
#####################

# Process CTRL
con = data.frame(Gene_symbol='CTRL', countsFile[grep('^CTRL0', rownames(countsFile)),])
countsFile = countsFile[-c(grep('^CTRL0', rownames(countsFile))),]

# Process Intergenic
intergenic = data.frame(Gene_symbol='INTERGENIC', countsFile[grep('INTERGENIC', rownames(countsFile)),])
countsFile = countsFile[-c(grep('INTERGENIC', rownames(countsFile))),]

# Process Genes
countsFile = data.frame(Gene_symbol = rownames(countsFile), countsFile)
countsFile$Gene_symbol = substring(countsFile$Gene_symbol, 3) # Remove sg
countsFile$Gene_symbol = data.frame(do.call(rbind, strsplit(as.vector(countsFile$Gene_symbol), split = "_")))[,1]

# Add all data together
countsFile = rbind(countsFile, con, intergenic)

# Write mageck input file
dir.create('mageck_run')
mageck = data.frame(sgRNA = rownames(countsFile), countsFile)
colnames(mageck)[2] = 'gene'
write.table(mageck, 
            paste(paste(getwd(), "/mageck_run/mageck_input.crisprProcess", sep=''),method,format(Sys.time(), "%Y-%m-%d"), "txt", sep = "."), na="", 
            col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

######################
# GRAPHICS - Open document for writing figures
#####################

pdf(paste("crisprProcess_report", format(Sys.time(), "%Y-%m-%d"), "pdf", sep = "."), width=11, height=8.5)

######################
# GRAPHICS - Plot raw counts input
#####################

par(mar=c(12,5,4,2))
par(las=2)
barplot(colSums(countsFile[,-1]), 
	col=rainbow(ncol(countsFile[,-1])),
	main="Initial Library Sizes",
	cex.names=0.75)
	
###############################
# 4. Remove guide RNAs with < 100 reads in input and less than 3 sgRNAs
###############################

guidesToRemove = nrow(countsFile[countsFile[,inputColumn] < 100,])

countsFile = countsFile[countsFile[,inputColumn] >=100,]
genesToRemove = table(countsFile$Gene_symbol)
genesToRemove = names(genesToRemove[genesToRemove <3])
genesToKeep = setdiff(countsFile$Gene_symbol,genesToRemove)
countsFile = countsFile[countsFile$Gene_symbol %in% genesToKeep,]

# Write out statistics of genes and guides removed from data
fileConn = file('genes_guides_removed.txt')
writeLines(c('Genes and guides removed',
             paste('Number of guides removed: ', guidesToRemove, sep=''),
             paste('Number of genes removed: ', length(genesToRemove), sep='')
             ), fileConn)
close(fileConn)

###############################
# 5. Generate log2 ratio difference between each sample and the input sample
###############################

countsFileMat = countsFile[,-1]
countsFileMat = countsFileMat + 1

# Calculate ratios based on input column name
for (i in 1:ncol(countsFileMat)){
		countsFileMat[,i] = countsFileMat[,i] / countsFileMat[,inputColumn]
}
countsFileMat = log2(countsFileMat)
countsFileMat = countsFileMat[,-which(colnames(countsFileMat) == inputColumn)]

######################
# GRAPHICS - Plot log2 ratio boxplots
#####################

par(mar=c(12,5,4,2))
par(las=2)
boxplot(countsFileMat, 
	col=rainbow(ncol(countsFileMat)),
	main="Distribution of Log2(late/input) Ratios",
	ylab = 'log2(late/input)',
	cex.names=0.75)

######################
# GRAPHICS - # View drop out line plots
#####################

par(mar=c(12,5,4,2))
par(las=1)
plot(sort(countsFileMat[,1]), type='l', col=rainbow(ncol(countsFileMat))[1], lwd=2, xlab='sgRNA Rank', ylab='Raw CRISPR Score (no centering)', cex.axis=1.25, cex.lab=1.25, 
	main = "Distribution of sgRNA Effects")
legend(x="bottomright", title="Cell Line",
   legend=c(colnames(countsFileMat)[1:ncol(countsFileMat)], 'Random distribution'),
   col=c(rainbow(ncol(countsFileMat))[1:ncol(countsFileMat)], 'black'), lwd=2, lty=c(rep(1,ncol(countsFileMat)),2), 
   pch=NA, bty="n", cex=0.7)
for (i in 2:ncol(countsFileMat)){
	lines(sort(countsFileMat[,i]), type='l', col=rainbow(ncol(countsFileMat))[i], lwd=2)
}
lines(sort(runif(nrow(countsFileMat), min=min(countsFileMat[,1]), max=max(countsFileMat[,1]))), type='l', col='black', lwd=1, lty=2)		

###############################
# 6. Normalize based on user input: Quantile normalization followed by median centering or zmad scaling
###############################

if (method == 'quantile'){
	
	suppressPackageStartupMessages(library(preprocessCore))
	library(preprocessCore)
	countsFileMat = normalize.quantiles(as.matrix(countsFileMat),copy=FALSE)
	myMedian = mean(apply(countsFileMat, 2, median))
	countsFileMat = countsFileMat - myMedian

} else if (method == 'zmad'){

	for (i in 1:ncol(countsFileMat)){
	countsFileMat[,i] = (countsFileMat[,i] - median(countsFileMat[,i])) / mad(countsFileMat[,i])
	}

} else {
	stop('Normalization method has not been selected: quantile or zmad only')
}

######################
# GRAPHICS - # View of sample distributions after normalization
#####################

par(mar=c(15,4,4,2))
par(las=2)
boxplot(countsFileMat,
	ylab="Normalized Score",
	col=rainbow(ncol(countsFileMat)),
	main=paste("Library representation after ", method, " normalization", sep=""))

###############################
# 7. Output normalized sgRNA data to text file 
###############################

output = data.frame(Guide_name = rownames(countsFileMat), Gene_symbol = countsFile$Gene_symbol, countsFileMat)
write.table(output, 
		paste("guide_level.crisprProcess",method,format(Sys.time(), "%Y-%m-%d"), "txt", sep = "."), na="", 
		col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

###############################
# 8. Calculate mean per gene data of normalized sgRNA values
###############################

countsFileMatMean = data.frame(Gene_symbol = countsFile$Gene_symbol, countsFileMat)
countsFileMatMean = aggregate(countsFileMatMean, by=list(countsFileMatMean$Gene_symbol), mean)
countsFileMatMean = countsFileMatMean[,-2]
colnames(countsFileMatMean)[1] = 'Gene_symbol'
write.table(countsFileMatMean, 
		paste("gene_level.crisprProcess",method,format(Sys.time(), "%Y-%m-%d"), "txt", sep = "."), na="", 
		col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

###############################
# 9. Write out gene symbols in the top 5% and bottom 5% of all samples
###############################

# Bottom 5%
fivePercent = countsFileMatMean
rownames(fivePercent) = fivePercent[,1]
fivePercent = fivePercent[,-1]
fivePercentGenes = as.character()

for (i in 1:ncol(fivePercent)){
  df = fivePercent[,i, drop=F]
  fivePercentGenes = c(fivePercentGenes,rownames(df[df[,1] <= quantile(df[,1], probs = 0.05),,drop=F]))
}

bottomFive = as.data.frame(table(fivePercentGenes))
bottomFive = bottomFive[order(bottomFive[,2], decreasing = TRUE),]
colnames(bottomFive) = c('Gene_symbol', 'Bottom_Five_Percentile_Count')
write.table(bottomFive, 
            paste("bottom_five_percent_genes",method,format(Sys.time(), "%Y-%m-%d"), "txt", sep = "."), na="", 
            col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


######################
# GRAPHICS - # Close graphics plotting
#####################

dev.off()
	
###############################
# 10. Generate shell script to run mageck for all comparisons
###############################		

mageckNames = colnames(mageck[,-which(colnames(mageck) %in% c('sgRNA', 'gene', inputColumn))])
mageckCommands = as.character()

# Generate list to store
for (i in 1:length(mageckNames)){
  mageckCommands = c(mageckCommands,paste('mageck test -k ', paste("mageck_input.crisprProcess",method,format(Sys.time(), "%Y-%m-%d"), "txt", sep = "."),' -t ',mageckNames[i],' -c ',inputColumn,' -n ',mageckNames[i], sep=''))
}

# Write shell script commands
fileConn = file(paste(getwd(), '/mageck_run/mageck_run.sh', sep=''))
writeLines(c('module load python/2.7.8',
             mageckCommands), fileConn)
close(fileConn)

print("CRISPR processing is complete.")


