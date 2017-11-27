#!/usr/bin/env python
#
# slice_alignment.py created 2017-11-26

'''slice_alignment.py  last modified 2017-11-26

slice_alignment.py -a matrix.phy -i 135,289,455 -o subset_matrix.phy

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm
'''

import sys
import os
import argparse
import time
import gzip
from collections import Counter
from Bio import AlignIO

def get_indices(indexfile):
	'''read comma-delimited or line separated file and return a dict where key is position'''
	indexdict = {} # key is int of index, value is True
	for line in open(indexfile,'r'):
		line = line.strip()
		if line:
			positions = line.split(",") # split "136,299,468,486"]
			for pos in positions: # keep positions as native numbers, 1-100, not 0-99
				indexdict[int(pos)] = True
	print >> sys.stderr, "# read {} positions from {}".format(len(indexdict), indexfile), time.asctime()
	return indexdict

def string_to_indices(indexstring):
	indexdict = {} # key is int of index, value is True
	positions = line.split(",") # split "136,299,468,486"]
	for pos in positions: # keep positions as native numbers, 1-100, not 0-99
		indexdict[int(pos)] = True
	print >> sys.stderr, "# read {} positions from input".format( len(indexdict) ), time.asctime()
	return indexdict

def slice_alignment(fullalignment, alignformat, indexdict):
	'''read alignment, and return a new alignment where partitions are reordered from highest coverage to lowest'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()

	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	num_species = len(alignedseqs)
	print >> sys.stderr, "# alignment has {} taxa with {} positions".format(num_species, alignedseqs.get_alignment_length()), time.asctime()

	newalign = alignedseqs[:,0:0] # start with blank alignment

	print >> sys.stderr, "# slicing alignment", time.asctime()
	for pos in range(alignedseqs.get_alignment_length()):
		if pos+1 in indexdict:
			alignpart = alignedseqs[:, pos:pos+1 ] # alignment of each partition only
			#print >> sys.stderr, "keeping pos {}".format(pos+1)
			newalign += alignpart
	print >> sys.stderr, "# alignment has {} taxa with {} sites".format( len(newalign), newalign.get_alignment_length() ), time.asctime()
	return newalign

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-i','--index', help="comma-separated list of positions (as numbers, no spaces) or file with numbers on each line")
	parser.add_argument('-o','--output', help="output name for new alignment", required=True)
	args = parser.parse_args(argv)

	if os.path.isfile(args.index):
		indexdict = get_indices(args.index)
	else: # assume string
		indexdict = string_to_indices(args.index)
	subalignment = slice_alignment(args.alignment, args.format, indexdict)
	AlignIO.write(subalignment, args.output, args.format)
	print >> sys.stderr, "# Subalignment written to {}".format(args.output), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
