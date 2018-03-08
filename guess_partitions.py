#!/usr/bin/env python
#
# guess_partitions.py  created 2018-01-13

'''guess_partitions.py  last modified 2018-03-07
  check a supermatrix and guess partitions based on changes of occupancy

guess_partitions.py -a supermatrix.phy -f phylip-relaxed

  set minimum number of changes to detect the site using -p
    increase the number until a plausible number of partitions becomes evident

  if missing genes have a unique character, like ?, use -q or -x
    as this is substantially more reliable than regular gaps

  output consists of three columns, of partition, length, and number of gaps
    between that site and the site immediately following
    so a line like:

1:360   360     25

  means the partition is predicted from site 1 to 360, where 25 letters or
    gaps change between site 360 and 361

  some manual editing may be required due to actual gaps

  partitions can then be easily extracted with

cut -f 1 matrix_pre-partitions.tab > partitions.txt
'''

import sys
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO

def is_gap(site, gapset):
	if site in gapset:
		return 1
	else:
		return 0

def check_alignments(fullalignment, alignformat, gapset, minlength=5, peakthreshold=2):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gaphisto = defaultdict(int)
	sitecov = {} # key is site number, value is coverage

	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	allength = alignedseqs.get_alignment_length()
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( numtaxa, allength )

	switchbysite = {} # key is site number, values is number of switches
	peakcount = 0
	lastpos = 1

	for i in range(allength-1): # i is python index
		changecount = 0
		firstsite = alignedseqs[:,i]
		secondsite = alignedseqs[:,i+1]
		for a1,a2 in zip(firstsite, secondsite):
			if is_gap(a1, gapset) + is_gap(a2, gapset) == 1:
				changecount += 1
		switchbysite[i] = changecount
		if changecount > peakthreshold:
			# i + 2 (add 1 for last index, add 1 to include the site itself) - lastpos
			partlength = i+2-lastpos
			if partlength < minlength: # if shorter than 5, skip and merge with next
				continue
			peakcount += 1
			print >> sys.stdout, "{}:{}\t{}\t{}".format(lastpos, i+1, partlength, changecount)
			lastpos = i+2 # offset by two, one for python index, one for next site
	else: # print last partition
		partlength = i+1-lastpos
		print >> sys.stdout, "{}:{}\t{}\t{}".format(lastpos, i+1, partlength, changecount)
	print >> sys.stderr, "# Found {} peaks, max change count was {}".format( peakcount, max(switchbysite.values()) )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-p','--peak-minimum', default=2, type=int, help="minimum peak height [2]")
	parser.add_argument('--len-minimum', default=5, type=int, help="minimum length to merge with next partition [5]")
	parser.add_argument('-q','--question-missing', action="store_true", help="use only question mark ? as missing gene, incompatible with -x")
	parser.add_argument('-x','--x-missing', action="store_true", help="use only X as missing gene, incompatible with -q")
	args = parser.parse_args(argv)

	GAPSET = "-?X"
	if args.question_missing: # reassign as only question mark
		GAPSET = "?"
	elif args.x_missing:
		GAPSET = "X"

	check_alignments(args.alignment, args.format, GAPSET, args.len_minimum, args.peak_minimum)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
