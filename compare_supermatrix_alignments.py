#!/usr/bin/env python
#
# compare_supermatrix_alignments.py created 2017-11-08

'''compare_supermatrix_alignments.py v1.0 2017-11-09
tool to compare different outputs from add_taxa_to_align.py

compare_supermatrix_alignments.py -1 matrix1.aln -2 matrix2.aln -p partitions.txt > comparison_matrix.tab

    for alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    partitions and all taxa must be the same for both alignments, meaning only
    differ by gene presence/absence, paralogs, or splice variation

    large matrices can be gzipped, as .gz

    comparison matrix file is printed to stdout, where matrix has five values for:
    same in both (2), different (1), and absent (0), found in #1 (-1), found in #2 (-2)

    this matrix can be plotted using the R script, where the output
    PDF is automatically renamed
Rscript draw_matrix_occupancy.R comparison_matrix.tab
'''

import sys
import argparse
import time
import gzip
from itertools import izip
from collections import defaultdict,Counter
from Bio import AlignIO

def get_partitions(partitionfile):
	'''read comma-delimited partition information and return a list of tuples'''
	partitions = [] # list of tuples of intervals
	for line in open(partitionfile,'r'):
		line = line.strip()
		if line:
			blocks = line.split(",") # split "1:136,137:301,..." into ['1:136', '137:301',...]
			for block in blocks:
				alignindex = tuple( int(i) for i in block.split(":") ) # split '1:136' into ( 1,136 )
				partitions.append(alignindex)
	print >> sys.stderr, "# read {} partitions from {}".format(len(partitions), partitionfile), time.asctime()
	return partitions

def check_alignments(alignfile1, alignfile2, format1, format2, partitions):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	comparisoncounts = defaultdict(int)
	occmatrix = []

	if alignfile1.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(alignfile1), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(alignfile1), time.asctime()
	alignedseqs1 = AlignIO.read(opentype(alignfile1), format1)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( len(alignedseqs1), alignedseqs1.get_alignment_length() )

	if alignfile2.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(alignfile2), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(alignfile2), time.asctime()
	alignedseqs2 = AlignIO.read(opentype(alignfile2), format2)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( len(alignedseqs2), alignedseqs2.get_alignment_length() )

	for seqrec in alignedseqs1: # all species must be counted once
		occmatrix.append([seqrec.id]) # initate each species as a new list with ID

	descriptiondict = { 0:"absent in both", 1:"different between the two", 2:"same in both", -1:"only in alignment 1", -2:"only in alignment 2"}
	for part in partitions:
		alignpart1 = alignedseqs1[:, part[0]-1:part[1] ] # alignment of each partition only
		alignpart2 = alignedseqs2[:, part[0]-1:part[1] ] # alignment of each partition only
		for i, (seqrec1, seqrec2) in enumerate( izip(alignpart1,alignpart2) ):
			species = seqrec1.id
			if species != seqrec2.id:
				print >> sys.stderr, "# WARNING: {} IS NOT SAME SPECIES AS {}".format( seqrec1.id, seqrec2.id )
			seqlen1 = len(seqrec1.seq)
			seqlen2 = len(seqrec2.seq)
			lettercounts1 = Counter(str(seqrec1.seq).replace("X","-").replace("?","-"))
			lettercounts2 = Counter(str(seqrec2.seq).replace("X","-").replace("?","-"))
			occupancyscore = 0 # by default is absent in both, reassign in all other cases
			if lettercounts1["-"] == seqlen1:
				if lettercounts2["-"] == seqlen2:
					occupancyscore = 0 # both all gaps, keep as 0
				else: # seq 1 is all gaps, seq 2 is present
					occupancyscore = -2
			else: # seq 1 is present
				if lettercounts2["-"] == seqlen2: # seq 2 is absent
					occupancyscore = -1
				else: # seq 2 is not all gaps
					if str(seqrec1.seq) == str(seqrec2.seq): # seqs are same
						occupancyscore = 2
					else: # seqs are different somehow
						occupancyscore = 1
			comparisoncounts[occupancyscore] += 1
			occmatrix[i].append(str(occupancyscore))
	for k,v in comparisoncounts.iteritems():
		print >> sys.stderr, "# {} sequences with score {}, meaning {}".format(v,k, descriptiondict[k])
	return occmatrix

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-1','--alignment1', help="supermatrix alignment 1")
	parser.add_argument('-2','--alignment2', help="supermatrix alignment 2, for comparison")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-D','--matrix-delimiter', default="\t", help="delimiter for matrix file [default is tab]")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	occmatrix = check_alignments(args.alignment1, args.alignment2, args.format, args.format, partitions)
	numparts = len(partitions)

	if occmatrix:
		print >> sys.stderr, "# writing matrix", time.asctime()
		# generate header line
		headerline = ["Species"] + ["{}-{}".format(*part) for part in partitions]
		print >> wayout, args.matrix_delimiter.join(headerline)
		# print occupancy by each species
		for occbysplist in occmatrix:
			print >> wayout, args.matrix_delimiter.join(occbysplist)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
