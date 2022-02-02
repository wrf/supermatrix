#!/usr/bin/env python
#
# trim_alignment_by_coverage.py  created 2017-11-27

'''trim_alignment_by_coverage.py  last modified 2022-02-02
    remove all gaps or sites with only a few taxa (default 1, set by -c)

trim_alignment_by_coverage.py -a alignment.phy -f phylip-relaxed

    can work with multiple files:

trim_alignment_by_coverage.py -a alignments/*.aln > trim_stats.tab

    output files are automatically renamed with .trim
    so alignment.phy would become alignment.phy.c1trim
'''

import sys
import argparse
import time
import gzip
from collections import Counter
from Bio import AlignIO

def trim_alignment(fullalignment, alignformat, covcutoff):
	'''read single alignment, print out new alignment with low-coverage sites removed'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write( "# reading alignment {} as gzipped  {}\n".format(fullalignment, time.asctime() ) )
	else: # otherwise assume normal open
		opentype = open
		sys.stderr.write( "# reading alignment {}  {}\n".format(fullalignment, time.asctime() ) )
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	allength = alignedseqs.get_alignment_length()
	sys.stderr.write( "# Alignment contains {} taxa for {} sites, including gaps\n".format( numtaxa, allength ) )

	allgaps = 0
	trimmedsites = 0
	keptsites = 0
	newalign = alignedseqs[:,0:0] # start with blank alignment
	for i in range(allength):
		aacounter = Counter( alignedseqs[:,i] ) # count all characters including gaps
		gapchars = aacounter["-"] + aacounter["X"] + aacounter["?"] # sum of gaps
		coverage = numtaxa - gapchars
		if coverage > covcutoff:
			keptsites += 1
			newalign += alignedseqs[:,i:i+1]
		else:
			if gapchars==numtaxa:
				allgaps += 1
			else:
				trimmedsites += 1
	newalignname = "{}.c{}trim".format(fullalignment, covcutoff)
	print("# {}\t{} gaps\t{} trimmed\t{} kept".format( newalignname, allgaps, trimmedsites, keptsites), file=sys.stdout)
	AlignIO.write(newalign, newalignname, alignformat)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="untrimmed alignment (can be many as *.aln)")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-c','--cutoff', type=int, default=1, help="coverage cutoff [1]")
	args = parser.parse_args(argv)

	for alignment in args.alignments:
		trim_alignment(alignment, args.format, args.cutoff)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
