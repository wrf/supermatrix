#!/usr/bin/env python
#
# reorder_alignment_by_tree.py created 2017-11-03

'''reorder_alignment_by_tree.py  last modified 2021-11-10

reorder_alignment_by_tree.py -a matrix.phy -T tree.nex -f phylip-relaxed > reordered.aln

    large matrices can be gzipped, as .gz

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm
'''

import sys
import os
import argparse
import time
import gzip
from collections import defaultdict,Counter
from Bio import AlignIO
from Bio import Phylo

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-T','--matrix-tree', help="optional Nexus-format tree to reorder matrix")
	args = parser.parse_args(argv)

	if args.alignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write( "# reading alignment {} as gzipped  {}\n".format(args.alignment, time.asctime() ) )
	else: # otherwise assume normal open
		opentype = open
		sys.stderr.write( "# reading alignment {}  {}\n".format(args.alignment, time.asctime() ) )
	alignedseqs = AlignIO.read(opentype(args.alignment), args.format)
	sys.stderr.write( "# Alignment contains {} taxa for {} sites, including gaps\n".format( len(alignedseqs), alignedseqs.get_alignment_length() ) )

	seqdict = {}
	for seqrec in alignedseqs:
		seqdict[seqrec.id] = seqrec

	sys.stderr.write( "# reading tree {}  {}\n".format(args.matrix_tree, time.asctime() ) )
	tree = Phylo.read(args.matrix_tree,"nexus")
	cladecount = 0
	for clade in tree.get_terminals():
		cleanname = str(clade.name).replace("'","").replace('"','')
		wayout.write( seqdict[cleanname].format("fasta") )
		cladecount += 1
	sys.stderr.write( "# wrote {} taxa  {}\n".format(cladecount, time.asctime() ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
