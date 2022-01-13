#! /usr/bin/env python
# join_alignments.py v1.0 created 2017-08-15

'''
join_alignments.py v1.1  last modified 2022-01-13

    specify multiple alignments with -a
    species names are expected to be the same, perhaps separated by an
      optional delimiter with -d
    for instance, Mnemiopsis_leidyi@ML172115a
      can be split with -d "@"

join_alignments.py -a hehenberger2017_alignments/* -d "@" -u hehenberger2017_supermatrix.fasta

    joined alignment -u will always be written in fasta format

    if rejoining alignments from add_taxa_to_align, after manual editing
    use -A to keep partition order
'''

import sys
import argparse
import time
import re
from Bio import AlignIO
from collections import OrderedDict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', nargs="*", help="alignment files", required=True)
	parser.add_argument('-A','--add-taxa-order', action="store_true", help="reorder partitions based on output names from add_taxa_to_align.py")
	parser.add_argument('-d','--delimiter', help="optional delimiter for sequence names")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-s','--sort', action="store_true", help="sort sequences alphabetically for output")
	parser.add_argument('-u','--supermatrix', help="name for supermatrix output", required=True)
	args = parser.parse_args(argv)

	runningsum = 0
	aligncounter = 0
	superprotsbytaxa = OrderedDict() # keys are taxa, values are concatenated proteins
	partitionlist = []

	problemtaxa = 0
	problemalignments = {}

	alignlist = args.alignments
	if args.add_taxa_order: # sort by partition number, assuming from add_taxa_to_align
		alignlist = sorted(alignlist, key=lambda x: int(re.search("(\d+)_(\d+)_part.aln",x).group(1)))

	for alignfile in alignlist:
		aligncounter += 1
		alignment = AlignIO.read(alignfile, args.format)
		al_length = alignment.get_alignment_length()
		targetlength = runningsum + al_length
		print( "# {} contains {} taxa for {} sites".format( alignfile, len(alignment), al_length ), file=sys.stderr) 

		# reset current keys each new alignment
		# must be done to catch for taxa that were in earlier alignments but not in current one
		existingkeys = dict([ (k,True) for k in superprotsbytaxa.keys() ])

		for seqrec in alignment:
			if args.delimiter:
				taxon_id = str(seqrec.id).split(args.delimiter)[0]
			else:
				taxon_id = str(seqrec.id)
			if taxon_id in superprotsbytaxa:
				try:
					existingkeys.pop(taxon_id)
					superprotsbytaxa[taxon_id] += str(seqrec.seq)
				except KeyError:
					print(  "WARNING: {} OCCURS MORE THAN ONCE IN ALIGNMENT {}".format( taxon_id, alignfile ) , file=sys.stderr )
					problemtaxa += 1
					problemalignments[alignfile] = True
			else: # if that taxa is not yet in growing matrix, fill entry with gaps
				superprotsbytaxa[taxon_id] = "-" * runningsum
				superprotsbytaxa[taxon_id] += str(seqrec.seq)
		# if any existing taxa are not in the current alignment, place all gaps
		for notfoundkey in existingkeys:
			superprotsbytaxa[notfoundkey] += "-" * al_length
		# extend the current running sum of lengths and update partitions
		partitionlist.append("{}:{}".format(runningsum+1, targetlength) )
		runningsum += al_length
	print( "# Finished parsing {} alignments for {} total taxa".format( aligncounter, len(superprotsbytaxa) ) , file=sys.stderr )

	### BUILD SUPERMATRIX
	with open(args.supermatrix,'w') as sm:
		if args.sort: # sort alphabetically
			for taxon in sorted( superprotsbytaxa.keys() ):
				print( ">{}\n{}".format(taxon, superprotsbytaxa[taxon]) , file=sm )
		else: # sequence order is arbitrary, based on whatever the first alignment had
			for taxon,sequence in superprotsbytaxa.items():
				print( ">{}\n{}".format(taxon, sequence) , file=sm )
		print( "# Supermatrix written to {}  {}".format( args.supermatrix, time.asctime() ), file=sys.stderr )
		with open("{}.partition.txt".format(args.supermatrix),'w') as pf:
			print( ",".join(partitionlist), file=pf )

	if problemalignments:
		print( "WARNING: {} ALIGNMENTS HAVE {} REDUNDANT TAXA".format( len(problemalignments), problemtaxa ) , file=sys.stderr )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
