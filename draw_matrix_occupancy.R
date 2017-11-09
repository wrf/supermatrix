#!/usr/bin/env Rscript
# generate figure from phylogenetic occupancy matrix
# v1 created 2017-08-03

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/supermatrix/philippe2009_occupancy_matrix.tab"
outputfile = gsub("([\\w/]+)\\....","\\1.pdf",inputfile,perl=TRUE)

# read table assuming tab delimited
occmatrix = read.table(inputfile,header=TRUE, sep="\t",row.names=1)

# get dimensions
taxa = dim(occmatrix)[1]
genes = dim(occmatrix)[2]

# set colors for 0, 1 and 2, as white, orange, and blue
occupancycolors = c("#fec44f","#addd8e","#FFFFFF","#ea4f12", "#477df6")

# generate matrix and reverse row order to place 1 at the top
m1 = matrix(unlist(occmatrix), ncol=taxa, byrow=TRUE)
m2 = m1[,ncol(m1):1]


# #
# COLOR TAXA BASED ON OCCUPANCY
# #
### SET THESE CONSTANTS TO CHANGE THE OCCUPANCY CUTOFF ###
OCCTHRESBYTAXA = 1.0
BADTHRESBYTAXA = 0.5
OCCTHRESBYGENE = 0.95

# calculated highest possible value per row as 2x number of genes
maxgscore = 2*genes
maxtscore = 2*taxa
# get list of taxa that have the maximum possible score, meaning 100% occupancy
print(paste("# coloring taxa with occupancy >=", OCCTHRESBYTAXA, sep=" " ))
completetaxa = taxa-which(apply(occmatrix, 1, sum)>=maxgscore*OCCTHRESBYTAXA)+1
# generate a color vector of all black
taxacolors = rep("#000000",taxa)
# reassign 100% complete taxa as green
taxacolors[completetaxa] = "#0aab30"
# reassign incomplete taxa as purple
print(paste("# coloring taxa with occupancy <", BADTHRESBYTAXA, sep=" " ))
badtaxa = taxa-which(apply(occmatrix, 1, sum)<maxgscore*BADTHRESBYTAXA)+1
taxacolors[badtaxa] = "#cd13c2"

# repeat for gene partitions to teal
print(paste("# coloring genes with occupancy >=", OCCTHRESBYGENE, sep=" " ))
completegenes = which(apply(occmatrix, 2, sum)>=maxtscore*OCCTHRESBYGENE)
genecolors = rep("#000000",genes)
genecolors[completegenes] = "#06ac7f"

# #
# GENERATE PDF OUTPUT #
# #
pdf(outputfile, width=4+genes%/%18, height=5+taxa%/%10)
par(mar=c(6,12,4,1))
# heatmap with order of matrix where first entry is at bottom
#image(x=1:genes, y=1:taxa, z=m1, col=occupancycolors, xlab="Genes", ylab="", axes=FALSE, main=inputfile )

# inverted order, to resemble viewing the table directly
# this should be considered "normal"
image(x=1:genes, y=1:taxa, z=m2, zlim=c(-2,2), col=occupancycolors, xlab="", ylab="", axes=FALSE, main=inputfile )
mtext("Gene partitions", 1, line=4, cex=1.3)
mtext(rev(row.names(occmatrix)),side=2,at=c(1:taxa),las=1, cex=0.9, col=taxacolors)
mtext(sub("X","",colnames(occmatrix)),side=1,at=c(1:genes),las=2, cex=0.7, col=genecolors)

# turn on XPD to allow drawing outside of normal frame
par(xpd=TRUE)
# if negative values are present, assume comparison script, and draw different legend
if (min(m2)<0) {
legend(-35,-1, legend=c("In m2", "In m1", "Absent","Mismatch","Same"), pch=22, pt.bg=occupancycolors, cex=1.1, ncol=2)
} else {
legend(-24,-1, legend=c("Absent","Incomplete","Present"), pch=22, pt.bg=occupancycolors[3:5], cex=1.1)
}

dev.off()

#