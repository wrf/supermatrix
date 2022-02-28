# plot homolog output logs
#
# last modified 2022-02-28
#
# based on output files from makehomologs.py
# https://bitbucket.org/wrf/sequences/src/master/makehomologs.py
#
# R script can be run as:
# Rscript plot_homolog_output_logs.R something.mh.log fasta_clusters.something.txt
#
# PDF will be automatically generated

args = commandArgs(trailingOnly=TRUE)

input_logfile = args[1]
#input_logfile = "~/project/lactobacillus_genomics/lactobacillus_v3.2022-02-25-092353.mh.log"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",input_logfile,perl=TRUE)

# read log, and then parse subparts of the log file
print(paste("Reading log file", input_logfile, ", writing to", outputfile))
log_data = read.table(input_logfile, stringsAsFactors = FALSE,skip=2, header=FALSE, sep="\t")

clustcounts = log_data[log_data$V1=="N_seqs",]
taxacounts = log_data[log_data$V1=="N_taxa",]
count_by_sp = log_data[log_data$V1=="N_genes",]

input_cluster_file = args[2]
#input_cluster_file = "~/project/lactobacillus_genomics/fasta_clusters.lactobacillus_v3.txt"

print(paste("Reading orthogroups from", input_cluster_file, ", writing to", outputfile))
cluster_data = read.table(input_cluster_file, header=FALSE, sep="\t")

num_sequences = as.integer(cluster_data$V2)
num_taxa = as.integer(cluster_data$V3)
min_per_taxon = as.integer(cluster_data$V4)
max_per_taxon = as.integer(cluster_data$V6)


# make PDF
print("Making PDF")
pdf(file=outputfile, width=8, height=11)
# page 1
par(mfrow=c(2,1), mar=c(4.5,4.5,2.5,1))
# upper plot
plot(clustcounts$V2, clustcounts$V3, type="n", xlab="Cluster size", ylab="Number of clusters",
     main=paste("Sizes of",sum(as.integer(clustcounts$V3)),"candidate clusters"),
     cex.axis=1.4, cex.lab=1.4)
lines(clustcounts$V2, clustcounts$V3 , lwd=5, col="#016c5999" )
lines(taxacounts$V2, taxacounts$V3 , lwd=5, col="#fc4e2a99" )
legend( "topright", legend=c("N sequences", "N taxa"), 
        col=c("#016c59","#fc4e2a"), lwd=8, cex=1.3 )
# lower plot
species_count = nrow(count_by_sp)
bar_las = ifelse( species_count < 35, 3, 0)
if (species_count < 35){
  barplot_names = count_by_sp$V2
  par(mar=c(8,4.5,2.5,1))
} else {
  barplot_names = NULL
  par(mar=c(3,4.5,2.5,1))
}
barplot( as.integer(count_by_sp$V3) , ylab="Number of genes in clusters", 
         main=paste("Total counts of genes in orthogroups for",species_count,"species"),
         col="#9ecae1", cex.axis=1.4, cex.lab=1.4 , names.arg=barplot_names, las=bar_las)

# end page 1

# page 2
par(mar=c(4.5,4.5,2.5,1))
is_avg_one = (num_sequences==num_taxa)
# upper plot
# orthologs are blue, paralogs are orange
point_color_vec = ifelse(is_avg_one,"#0000b344","#ec701444")
average_seqs_per_sp = num_sequences/num_taxa
plot(num_sequences, average_seqs_per_sp, xlab="Total sequences in cluster", ylab="Average seqs per species",
     main=paste("Average number of seqs per taxon of",nrow(cluster_data),"orthogroups"),
     cex.axis=1.4, cex.lab=1.4, #xlim=c(0,30), ylim=c(0,6),
     cex=2, col=point_color_vec, pch=16)
abline(a=0, b=1/species_count , lwd=1, lty=2, col="#00000044" ) # limit of number of species
abline(a=0, b=1/min(num_taxa) , lwd=1, lty=2, col="#00000044" ) # limit of -s , minimum taxa per cluster
text( max(num_sequences)*0.08, max(average_seqs_per_sp), ":", cex=1.3)
text( max(num_sequences)*0.09, max(average_seqs_per_sp), table(is_avg_one)[["TRUE"]], adj=0, col="#0000b3")
text( max(num_sequences)*0.07, max(average_seqs_per_sp), table(is_avg_one)[["FALSE"]], adj=1, col="#ec7014")
#text(1,max(num_sequences/num_taxa), "-s" )

# lower plot
# change color of paralogs to pink
point_color_vec = ifelse(is_avg_one,"#0000b344","#dd349744")
plot(num_sequences, max_per_taxon, xlab="Total sequences in cluster", ylab="Most seqs in a species",
     main=paste("Most seqs from a single taxon for",nrow(cluster_data),"orthogroups"),
     cex.axis=1.4, cex.lab=1.4, #xlim=c(0,30), ylim=c(0,12),
     cex=2, col=point_color_vec, pch=16)
abline(a=1-min(num_taxa), b=1, lwd=1, lty=2, col="#00000044" ) # limit of -s
abline(v=min(num_sequences), lwd=1, lty=2, col="#00000044" ) # vertical line limit of -z
text( max(num_sequences)*0.08, max(max_per_taxon), ":", cex=1.3)
text( max(num_sequences)*0.09, max(max_per_taxon), table(is_avg_one)[["TRUE"]], adj=0, col="#0000b3")
text( max(num_sequences)*0.07, max(max_per_taxon), table(is_avg_one)[["FALSE"]], adj=1, col="#dd3497")

dev.off()

print("Done")


#