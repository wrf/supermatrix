#!/usr/bin/env python
#
# coverage_by_site.py  created 2017-11-27

'''coverage_by_site.py  last modified 2017-11-27

coverage_by_site.py -a supermatrix.phy
'''

import sys
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO

def check_alignments(fullalignment, alignformat, printbysite):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gaphisto = defaultdict(int)

	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	allength = alignedseqs.get_alignment_length()
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( len(alignedseqs), allength )

	for i in range(allength):
		aacounter = Counter( alignedseqs[:,i] ) # count all characters including gaps
		gapchars = aacounter["-"] + aacounter["X"] + aacounter["?"]
		gaphisto[ gapchars ] += 1

	print >> sys.stderr, "numGaps\tnumSites\tcoverage"
	for k in sorted(gaphisto.keys()):
		print >> sys.stderr, "{}\t{}\t{}".format( k, gaphisto[k], len(alignedseqs)-k )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-m','--matrix-out', help="name for optional matrix-occupancy output file")
	parser.add_argument('-D','--matrix-delimiter', default="\t", help="delimiter for matrix file [default is tab]")
	args = parser.parse_args(argv)

	check_alignments(args.alignment, args.format, args.matrix_out)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
