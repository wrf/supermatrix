#!/usr/bin/env Rscript
# generate barplot from phylogenetic occupancy matrix for multiple matrices
# v1 created 2018-03-05

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/supermatrix/matrix/philippe2009_percent_matrix.tab"
outputfile = gsub("([\\w/]+)\\....","\\1.barplot.pdf",inputfile,perl=TRUE)

# read table assuming tab delimited
occmatrix = read.table(inputfile,header=TRUE, sep="\t",row.names=1)

# get dimensions
taxa = dim(occmatrix)[1]
genes = dim(occmatrix)[2]

sumsbyspecies = rowSums(occmatrix)/genes/100
meanoccupancy = mean(sumsbyspecies)

if (genes < 2) {
xaxislabel = paste("Fraction of occupied sites")
} else {
xaxislabel = paste("Fraction of occupied partitions ( n =",genes,")")
}

# #
# COLOR TAXA BASED ON OCCUPANCY
# #
### SET THESE CONSTANTS TO CHANGE THE OCCUPANCY CUTOFF ###

occupancycolors = c("#FFFFFF",colorRampPalette(c("#ea4f12", "#477df6"))(100))

# #
# GENERATE PDF OUTPUT #
# #
pdf(outputfile, width=10, height=5+taxa%/%10)
par(mar=c(4.5,15,4,1.2))
# inverted order, to resemble viewing the table directly
# this should be considered "normal"
barplot(rev(sumsbyspecies), xlim=c(0,1), horiz=TRUE, axes=FALSE, main=inputfile, names=rev(row.names(occmatrix)), border=NA, col=occupancycolors[round(rev(sumsbyspecies),digits=2)*100+1], las=1 )
axis(1, cex.axis=1.4)
mtext( xaxislabel , 1, line=3, cex=1.3)
mtext( paste("n taxa =",taxa) , side=3, line=1, at=c(-0.1), cex=1.3 )
lines(c(meanoccupancy,meanoccupancy),c(1,taxa*1.2), lwd=3, lty=2)
dev.off()

#
