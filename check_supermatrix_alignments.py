#!/usr/bin/env python
#
# check_supermatrix_alignments.py created 2017-03-13

'''check_supermatrix_alignments.py v1.0 2017-03-27
tool to quickly check for abnormal sequences in fasta alignments

checknogalignments.py -a matrix.phy -p partitions.txt

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    output consists of tab delimited fields:
Species  Partitions  Number-missing  Percent-missing  Number-partial  Percent-partial
'''

import sys
import os
import argparse
import time
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

def check_alignments(fullalignment, alignformat, partitions):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gapdict = {} # must set all values to zero in order to not skip full taxa
	halfgapdict = defaultdict(int)
	alignedseqs = AlignIO.read(fullalignment, alignformat)
	for seqrec in alignedseqs: # all species must be counted once
		gapdict[seqrec.id] = 0
	for part in partitions:
		alignpart = alignedseqs[:, part[0]-1:part[1] ]
		for seqrec in alignpart:
			species = seqrec.id
			seqlen = len(seqrec.seq)
			lettercounts = Counter(str(seqrec.seq).replace("X","-"))
			if lettercounts["-"] == seqlen: # seq is all gaps, so no seq
				gapdict[species] += 1
			elif lettercounts["-"] >= seqlen/2:
				halfgapdict[species] += 1
	print >> sys.stderr, "# split {} taxa alignment by partitions".format(len(gapdict)), time.asctime()
	return gapdict, halfgapdict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	gapdict, halfgaps = check_alignments(args.alignment, args.format, partitions)
	numparts = len(partitions)

	if args.header:
		#                  0           1           2       3
		print >> wayout, "Species\tPartitions\tMissing\tM%\tPartial\tP%"

	for k,v in sorted(gapdict.iteritems(), reverse=True, key=lambda x: x[1]):
		print >> wayout, "{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}".format(k, numparts, v, v*100.0/numparts, halfgaps[k], halfgaps[k]*100.0/numparts)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
