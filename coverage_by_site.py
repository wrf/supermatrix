#!/usr/bin/env python
#
# coverage_by_site.py  created 2017-11-27
# python3 update  2022-04-26

'''coverage_by_site.py  last modified 2022-04-26
  quick diagnostic of coverage
  reports histogram of number of sites with C coverage to stdout

coverage_by_site.py -a supermatrix.phy -f phylip-relaxed

  fields are:
num-missing-taxa    count-of-sites    num-taxa-with-site
    num-missing-taxa = M, number of gaps in that column
    count-of-sites = number of sites with missing-taxa M and coverage C
    num-taxa-with-site = coverage, by N-taxa - num-missing-taxa

  optional matrix of coverage by site can be given with -m

  sites are given as numbers (starting with 1)

coverage_by_site.py -a supermatrix.aln -m matrix_cov_by_site.tab
'''

import sys
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO

def check_alignments(fullalignment, alignformat, makeheader, trimlength):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	gaphisto = defaultdict(int)
	sitecov = {} # key is site number, value is coverage

	if fullalignment.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# reading alignment {} as gzipped  {}\n".format(fullalignment, time.asctime() ) )
	else: # otherwise assume normal open
		opentype = open
		sys.stderr.write("# reading alignment {}  {}\n".format(fullalignment, time.asctime() ) )

	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	numtaxa = len(alignedseqs)
	allength = alignedseqs.get_alignment_length()
	sys.stderr.write("# Alignment contains {} taxa for {} sites, including gaps\n".format( numtaxa, allength ) )

	constcounter = 0
	gapsum = 0
	charsum = 0

	effective_sites = allength
	if trimlength:
		effective_sites = trimlength
		sys.stderr.write( "# Using only first {} sites\n".format( trimlength ) )

	gap_only_sites = []

	for i in range(effective_sites):
		aacounter = Counter( alignedseqs[:,i] ) # count all characters including gaps
		gapchars = aacounter["-"] + aacounter["X"] + aacounter["?"] # sum of gaps
		gaphisto[ gapchars ] += 1 # key is sum of gaps, value is number of sites with n gaps
		gapsum += gapchars
		charsum += sum(aacounter.values())
		if gapchars==numtaxa:
			gap_only_sites.append(str(i+1))

		aaset = set(aacounter.keys())
		aasetnogaps = set(aaset)
		aasetnogaps.discard("-")

		numAAs = len(aasetnogaps)
		aasetstring = ",".join(["{}:{}".format(x,aacounter[x]) for x in list(aasetnogaps)])

		if aasetnogaps and len(aasetnogaps) == 1:
			constcounter += 1

		sitecov[i] = [ numtaxa-gapchars, numAAs, aasetstring ]

	if gap_only_sites:
		sys.stderr.write( "# WARNING: found {} sites containing only gaps\n{}\n".format(len(gap_only_sites), ", ".join(gap_only_sites) ) )

	if makeheader:
		sys.stdout.write("numGaps\tnumSites\tcoverage\n")
	for k in sorted(gaphisto.keys()):
		sys.stdout.write("{}\t{}\t{}\n".format( k, gaphisto[k], numtaxa-k ) )

	sys.stderr.write( "# {} total gaps out of {} characters ({:.2f}% total)\n".format( gapsum, charsum, 100.0 * gapsum/charsum) )
	averagegaps = 1.0*gapsum/effective_sites
	averagegaps_pct = 100*averagegaps/numtaxa
	sys.stderr.write( "# {:.2f} ({:.2f}% of {} taxa) average gaps per site, for {} sites\n".format( averagegaps, averagegaps_pct, numtaxa, effective_sites) )
	sys.stderr.write( "# {} constant sites (excluding any gaps at that site)\n".format(constcounter) )

	return sitecov

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-H','--header', action="store_true", help="include header line")
	parser.add_argument('-m','--matrix-out', help="name for optional matrix-occupancy output file")
	parser.add_argument('-D','--matrix-delimiter', default="\t", help="delimiter for matrix file [default is tab]")
	parser.add_argument('-t', '--trim', type=int, help="take only the first N letters")
	parser.add_argument('--aa-stats', action="store_true", help="print detailed output in matrix")
	args = parser.parse_args(argv)

	sitecovdict = check_alignments(args.alignment, args.format, args.header, args.trim)

	if args.matrix_out:
		with open(args.matrix_out,'w') as mo:
			for k in sorted(sitecovdict.keys()):
				if args.aa_stats:
					mo.write("{}{}{}\n".format( k+1, args.matrix_delimiter, args.matrix_delimiter.join( [str(x) for x in sitecovdict[k] ] ) ) )
				else:
					mo.write("{}{}{}\n".format( k+1, args.matrix_delimiter, sitecovdict[k][0] ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
