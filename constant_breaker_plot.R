#!/usr/bin/env Rscript
# generate figure from phylogenetic occupancy matrix
# v1 created 2018-04-16

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/est/supermatrix/v11/v11_linsi_c10trim_const_break_matrix.tab"
outputfile = gsub("([\\w/]+)\\....","\\1.pdf",inputfile,perl=TRUE)

# read table assuming tab delimited
occmatrix = read.table(inputfile,header=TRUE, sep="\t",row.names=1)

# get dimensions
taxa = dim(occmatrix)[1]
genes = dim(occmatrix)[2]

# generate matrix and reverse row order to place 1 at the top
m1 = matrix(unlist(occmatrix), ncol=taxa, byrow=TRUE)
m2 = m1[,ncol(m1):1]

#colmeds = apply(m2,1,median)
#medmatrix = matrix(rep(colmeds,taxa), ncol=taxa, byrow=TRUE)

#diffmatrix = m2-medmatrix

ERRTHRESBYTAXA = 1.0
BADTHRESBYTAXA = 0.33

occupancycolors = c("#F6F6F6",colorRampPalette(c("#a9d093","#5e3c99"))(11))

# generate a color vector of all black
taxacolors = rep("#000000",taxa)
# reassign incomplete taxa as purple
print(paste("# coloring taxa with occupancy <", BADTHRESBYTAXA, sep=" " ))
badtaxa = which(apply(m2, 2, function (x) sum(x[x == -1])) < -BADTHRESBYTAXA*genes)
taxacolors[badtaxa] = "#cd13c2"
print(paste("# coloring taxa with total errors >", genes, sep=" " ))
errorprone = which(apply(m2, 2, function (x) sum(x[x >= 0])) > ERRTHRESBYTAXA*genes)
taxacolors[errorprone] = "#ea4f12"

# #
# GENERATE PDF OUTPUT #
# #

pdf(outputfile, width=4+genes%/%18, height=5+taxa%/%10)
par(mar=c(6,12,4,1))
# inverted order, to resemble viewing the table directly
# this should be considered "normal"
image(x=1:genes, y=1:taxa, z=m2, zlim=c(-1,10), col=occupancycolors, xlab="", ylab="", axes=FALSE, main=inputfile )
mtext( paste("Gene partitions ( n =",genes,")") , 1, line=4, cex=1.3)
mtext( paste("n taxa =",taxa) , side=3, line=1, at=c(-10), cex=1.3 )
mtext(rev(row.names(occmatrix)),side=2,at=c(1:taxa),las=1, cex=0.9, col=taxacolors )
mtext(sub("X","",colnames(occmatrix)),side=1,at=c(1:genes),las=2, cex=0.7 )

# turn on XPD to allow drawing outside of normal frame
par(xpd=TRUE)
# if negative values are present, assume comparison script, and draw different legend
legend(-32,-1, legend=c("Absent", "0 score", "10 AAs" ), pch=22, pt.bg=occupancycolors[c(1,2,11)], cex=1.1)

dev.off()

#