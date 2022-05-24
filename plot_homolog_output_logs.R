# plot homolog output logs
#
# last modified 2022-05-23
#
# based on output files from makehomologs.py
# https://bitbucket.org/wrf/sequences/src/master/makehomologs.py
#
# R script can be run as:
# Rscript plot_homolog_output_logs.R something.mh.log fasta_clusters.something.tab
#
# PDF will be automatically generated

args = commandArgs(trailingOnly=TRUE)

input_logfile = args[1]
#input_logfile = "~/project/lactobacillus_genomics/lactobacillus_v3.2022-02-28-223648.mh.log"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",input_logfile,perl=TRUE)

# read log, and then parse subparts of the log file
print(paste("Reading log file", input_logfile, ", writing to", outputfile))
log_data = read.table(input_logfile, stringsAsFactors = FALSE,skip=2, header=FALSE, sep="\t")

clustcounts = log_data[log_data$V1=="N_seqs",]
taxacounts = log_data[log_data$V1=="N_taxa",]
genes_in_clusters_by_sp = log_data[log_data$V1=="N_genes",]
ortho_counts_by_sp = log_data[log_data$V1=="N_ortho",]


# these might be empty is blast processing step was skipped
total_prots_by_sp = log_data[log_data$V1=="T_genes",]
blast_hits_by_sp = log_data[log_data$V1=="N_bhits",]
prots_no_hits = as.integer(total_prots_by_sp$V3) - as.integer(blast_hits_by_sp$V3)
hits_no_clusters = as.integer(blast_hits_by_sp$V3) - as.integer(genes_in_clusters_by_sp$V3)
hits_w_clusters = as.integer(genes_in_clusters_by_sp$V3) - as.integer(ortho_counts_by_sp$V3)
ortho_clusters = as.integer(ortho_counts_by_sp$V3)

#
input_cluster_file = args[2]
#input_cluster_file = "~/project/lactobacillus_genomics/fasta_clusters.lactobacillus_v3.tab"

print(paste("Reading orthogroups from", input_cluster_file, ", writing to", outputfile))
# should autodetect if makehomologs.py used the -c option, making the header fasta_clusters.H.
cf_skip = grepl("fasta_clusters.H.", input_cluster_file)
# meaning skip the header line
cluster_data = read.table(input_cluster_file, header=cf_skip, sep="\t", stringsAsFactors = FALSE )


num_sequences = as.integer(cluster_data[,2])
num_taxa = as.integer(cluster_data[,3])
min_per_taxon = as.integer(cluster_data[,4])
med_per_taxon = as.integer(cluster_data[,5])
max_per_taxon = as.integer(cluster_data[,6])


# get vector for 1-1 orthologs
is_avg_one = (num_sequences==num_taxa)

ortho_only_sizes = table(num_sequences[is_avg_one])

# make PDF
print("Making PDF")
pdf(file=outputfile, width=8, height=11)
# page 1
par(mfrow=c(2,1), mar=c(4.5,4.5,2.5,1))
# upper plot
#first_xlim = min(rev(as.integer(clustcounts$V2))[cumsum(as.integer(rev(clustcounts$V3))) > 0.995*sum(as.integer(rev(clustcounts$V3)))])
plot(clustcounts$V2, clustcounts$V3, type="n", xlab="Cluster size", ylab="Number of clusters",
     main=paste("Sizes of",sum(as.integer(clustcounts$V3)),"candidate clusters"), # xlim=c(1,first_xlim),
     cex.axis=1.4, cex.lab=1.4 )
lines(clustcounts$V2, clustcounts$V3 , lwd=5, col="#016c5999" )
lines(taxacounts$V2, taxacounts$V3 , lwd=5, col="#fc4e2a99" )
lines( as.integer(unlist(dimnames(ortho_only_sizes))), ortho_only_sizes , lwd=5, col="#0000b399")
legend( "topright", legend=c("N sequences", "N taxa", "N orthologs"), 
        col=c("#016c59","#fc4e2a", "#0000b3"), lwd=8, cex=1.3 )
# lower plot
species_count = nrow(genes_in_clusters_by_sp)
bar_las = ifelse( species_count < 35, 3, 0)
if (species_count < 35){
  barplot_names = genes_in_clusters_by_sp$V2
  par(mar=c(8,4.5,2.5,1))
  bordcolor = "#000000"
} else {
  barplot_names = NULL
  par(mar=c(3,4.5,2.5,1))
  bordcolor = NA
}

if ( nrow(total_prots_by_sp)==nrow(genes_in_clusters_by_sp) ) {
  matrix_for_bars = t(as.matrix(data.frame( ortho_clusters, hits_w_clusters , hits_no_clusters , prots_no_hits )))
  barplot( matrix_for_bars , ylab="Number of genes in clusters", 
           main=paste("Total counts of prots, those with hits, those in orthogroups for",species_count,"species"),
           col=c("#0000b3", "#9ecae1","#8c6bb1","#7f2704"), border = bordcolor,
           cex.axis=1.4, cex.lab=1.4 , names.arg=barplot_names, las=bar_las)
  legend(x="bottom", legend=c("Single-copy", "Multi-copy", "Not clustered", "No blast"),
         col=c("#0000b3", "#9ecae1","#8c6bb1","#7f2704"), pch=15, pt.cex=2, horiz=TRUE, xpd = TRUE )
} else {
  barplot( as.integer(genes_in_clusters_by_sp$V3) , ylab="Number of genes in clusters", 
           main=paste("Total counts of genes in orthogroups for",species_count,"species"),
           col="#9ecae1", border = bordcolor,
           cex.axis=1.4, cex.lab=1.4 , names.arg=barplot_names, las=bar_las)
}
# end page 1


# page 2
par(mar=c(4.5,4.5,2.5,1))
# upper plot
point_color_vec = ifelse(is_avg_one,"#0000b344","#ed627844")
plot(num_sequences, num_taxa, xlab="Total sequences in cluster", ylab="Number of taxa in cluster",
     main=paste("Total sequences and taxon counts for",nrow(cluster_data),"orthogroups"),
     cex.axis=1.4, cex.lab=1.4, #xlim=c(0,30), ylim=c(0,12),
     cex=2, col=point_color_vec, pch=16)
abline(v=min(num_sequences), lwd=1, lty=2, col="#00000044" ) # limit of -z
text( max(num_sequences)*0.08, max(num_taxa), ":", cex=1.3)
text( max(num_sequences)*0.09, max(num_taxa), table(is_avg_one)[["TRUE"]], adj=0, col="#0000b3")
text( max(num_sequences)*0.07, max(num_taxa), table(is_avg_one)[["FALSE"]], adj=1, col="#ed6278")

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

# change color of paralogs to orange
has_more_than_two = (num_taxa > 2) # exclude clusters with only 2 taxa, as median is meaningless
median_vs_min = (med_per_taxon - min_per_taxon)[has_more_than_two]
max_vs_median = (max_per_taxon - med_per_taxon)[has_more_than_two]
point_color_vec = ifelse(is_avg_one,"#0000b344","#dd973444")[has_more_than_two]
plot( jitter(median_vs_min), jitter(max_vs_median), xlab="Median minus minimum", ylab="Maximum minus median",
     main=paste("Relative differences of min-med and med-max for",sum(has_more_than_two),"orthogroups"),
     cex.axis=1.4, cex.lab=1.4, #xlim=c(0,30), ylim=c(0,12),
     cex=2, col=point_color_vec, pch=16)

#d = data.frame(num_sequences, max_per_taxon)
#image(as.matrix(log10(table(d))), col=terrain.colors(10,rev = TRUE))

dev.off()

print("Done")


#