#!/usr/bin/env python
#
# check_supermatrix_alignments.py created 2017-03-13

'''check_supermatrix_alignments.py v1.1 2017-08-03
tool to quickly check for abnormal sequences in fasta alignments

checknogalignments.py -a matrix.phy -p partitions.txt

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    large matrices can be gzipped, as .gz

    output consists of tab delimited fields:
Species  Partitions  Number-missing  Percent-missing  Number-partial  Percent-partial

    for optional occupancy matrix file:
checknogalignments.py -a matrix.phy -p partitions.txt -m occupancy_matrix.tab
    where matrix consists of three values for:
    present (2), partial (1), and absent (0)
'''

import sys
import os
import argparse
import time
import gzip
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

def check_alignments(fullalignment, alignformat, partitions, makematrix=False):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gapdict = {} # must set all values to zero in order to not skip full taxa
	halfgapdict = defaultdict(int)

	occmatrix = [] if makematrix else None

	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open for GTF format
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	for seqrec in alignedseqs: # all species must be counted once
		gapdict[seqrec.id] = 0
		if occmatrix is not None:
			occmatrix.append([seqrec.id]) # initate each species as a new list with ID
	for part in partitions:
		alignpart = alignedseqs[:, part[0]-1:part[1] ] # alignment of each partition only
		for i,seqrec in enumerate(alignpart):
			species = seqrec.id
			seqlen = len(seqrec.seq)
			lettercounts = Counter(str(seqrec.seq).replace("X","-"))
			occupancyscore = 2 # by default is present, reassign if absent or partial
			if lettercounts["-"] == seqlen: # seq is all gaps, so no seq
				gapdict[species] += 1
				occupancyscore = 0 # set to 0 if all gaps
			elif lettercounts["-"] >= seqlen/2:
				halfgapdict[species] += 1
				occupancyscore = 1 # set to 1 if partial

			if occmatrix: # if building matrix, add that value to the matrix at sequence i
				occmatrix[i].append(str(occupancyscore))
	print >> sys.stderr, "# split {} taxa alignment by partitions".format(len(gapdict)), time.asctime()
	return gapdict, halfgapdict, occmatrix

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-m','--matrix-out', help="name for optional matrix-occupancy output file")
	parser.add_argument('-D','--matrix-delimiter', default="\t", help="delimiter for matrix file [default is tab]")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	gapdict, halfgaps, occmatrix = check_alignments(args.alignment, args.format, partitions, args.matrix_out)
	numparts = len(partitions)

	if args.matrix_out and occmatrix:
		print >> sys.stderr, "# writing matrix to {}".format(args.matrix_out), time.asctime()
		with open(args.matrix_out,'w') as mo:
			# generate header line
			headerline = ["Species"] + ["{}-{}".format(*part) for part in partitions]
			print >> mo, args.matrix_delimiter.join(headerline)
			# print occupancy by each species
			for occbysplist in occmatrix:
				print >> mo, args.matrix_delimiter.join(occbysplist)

	if args.header:
		#                  0           1           2       3
		print >> wayout, "Species\tPartitions\tMissing\tM%\tPartial\tP%"

	for k,v in sorted(gapdict.iteritems(), reverse=True, key=lambda x: x[1]):
		print >> wayout, "{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}".format(k, numparts, v, v*100.0/numparts, halfgaps[k], halfgaps[k]*100.0/numparts)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
