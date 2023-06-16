#!/usr/bin/env python
# make_ortholog_stats_table.py v1 2022-11-15

"""
make_ortholog_stats_table.py v1.0  last modified 2022-12-02

~/git/supermatrix/make_ortholog_stats_table.py -i clusters_lactobacillus_v4 > lactobacillus_cluster_stats_v4.tab

    or can glob with:

 ~/git/supermatrix/make_ortholog_stats_table.py -i clusters_lactobacillus_v4/* > lactobacillus_cluster_stats_v4.tab
"""

import sys
import os
import argparse
import time
from Bio import SeqIO
from glob import iglob

def parse_gene_names_from_db(proteinfile):
	"""from a fasta file, extract ID and gene name, description, or anything from the fasta header, 
	and return a dict, where key is ID and value is gene name, where available"""
	accession_to_gene = {}
	return accession_to_gene

def find_intermediate_ids(proteinfile):
	"""from a fasta file, get IDs and accessions, and return a dict of ID to accession"""
	prot_id_to_accession = {}
	return prot_id_to_accession

def get_gene_from_members(memberlist, accession_to_gene, prot_id_to_accession):
	"""for each gene, find the best or most common gene name"""
	usable_gene_names = {}
	best_name = ""
	return best_name

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', nargs="*", help="input files, or folder containing fasta files from makehomologs.py")
	parser.add_argument('-v','--verbose', action="store_true", help="display more output to stderr")
	parser.add_argument('--no-header', action="store_true", help="omit header line")
	args = parser.parse_args(argv)

	file_counter = 0
	protein_counter = 0
	duplicate_counter = 0

	accession_to_gene_dict = {}
	prot_id_to_accession_dict = {}

	if os.path.isdir(args.input[0]): #
		clust_list = os.path.join(args.input[0], "*") 
		clust_iter_obj = iglob(clust_list)
		sys.stderr.write("# Starting parsing on clusters in folder {}  {}\n".format( args.input[0], time.asctime() ) )
	elif os.path.isfile(args.input[0]): # first item is a file, hopefully the rest are
		clust_iter_obj = list(args.input)
		sys.stderr.write("# Starting parsing on clusters, first is {}  {}\n".format( args.input[0], time.asctime() ) )

	if not args.no_header:
		headerline = "cluster_ID\tgene_name\tnotes\tN_seqs\tmin_len\tmean_len\tmax_len\tmember_list\n"
		sys.stdout.write(headerline)

	for clustfile in clust_iter_obj:
		file_counter += 1
		# reset parameters for cluster
		best_gene_name = "NA"
		user_notes = ""
		seq_counter = 0
		min_length = 0
		max_length = 0
		member_list = {} # key is seqID, value is length
		for seqrec in SeqIO.parse(clustfile, "fasta"):
			protein_counter += 1
			seqlen = len(seqrec.seq)
			if seqrec.id not in member_list:
				member_list[seqrec.id] = seqlen
			else:
				duplicate_counter += 1
				sys.stderr.write( "WARNING: duplicate seq ID {} in cluster {}\n".format(seqrec.id, clustfile) )
		if member_list:
			min_length = min(member_list.values())
			max_length = max(member_list.values())
			mean_length = sum(member_list.values())/len(member_list.values())
			#best_gene_name = get_gene_from_members( list(member_list.keys()), accession_to_gene_dict, prot_id_to_accession_dict )

		outline = "{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\n".format(clustfile, best_gene_name, user_notes, len(member_list.values()), min_length, mean_length, max_length, " ; ".join( list(member_list.keys()) ) )
		sys.stdout.write(outline)

	sys.stderr.write("# Parsed {} clusters with {} proteins  {}\n".format( file_counter, protein_counter, time.asctime() ) )
	if duplicate_counter:
		sys.stderr.write("# Counted {} proteins with duplicate names\n".format( duplicate_counter ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

