#!/usr/bin/env python
#
# merge_taxa_pairs.py created 2018-01-16

'''merge_taxa_pairs.py v1.0 2018-01-16
    merge pairs of taxa keeping only complete proteins

merge_taxa_pairs.py -a matrix.phy -p partitions.txt -t Hsap,Hsap2 Mmus,Mmus2 > merged_taxa.aln

    merged sequences are written to stdout in fasta format
    this can be quickly merged with the original sequences

cat merged_taxa.aln matrix.aln > matrix_w_merged_taxa.aln

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    large matrices can be gzipped, as .gz
'''

import sys
import argparse
import time
import gzip
from itertools import izip
from collections import Counter
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

def merge_pairs(fullalignment, alignformat, partitions, taxonpairs, verbose):
	'''read large alignment, return the dict where key is species and value is number of gap-only sequences'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	print >> sys.stderr, "# Alignment contains {} taxa for {} sites, including gaps".format( len(alignedseqs), alignedseqs.get_alignment_length() )

	seqdict = {} # key is seqrec.id, value is seqrec
	for seqrec in alignedseqs: # need to make a dictionary to call species directly
		seqdict[seqrec.id] = seqrec

	for pair in taxonpairs:
		outcome_counter_dict = { "empty":0, "1st-yes-2nd-none":0, "1st-none-2nd-yes":0, "both-same":0, "1st-yes-2nd-partial":0, "1st-partial-2nd-yes":0, "both-different":0, "merged-partials":0}
		taxon1, taxon2 = pair.split(",")
		merge_id = taxon1 # by default, the merge uses name of taxon 1
		merge_string = ""
		print >> sys.stderr, "# merging {} and {}".format(taxon1, taxon2)
		for part in partitions:
			partlength = part[1] - part[0] + 1
			seq1part = str(seqdict[taxon1].seq)[part[0]-1:part[1]].replace("X","-").replace("?","-")
			seq2part = str(seqdict[taxon2].seq)[part[0]-1:part[1]].replace("X","-").replace("?","-")
			lettercounts1 = Counter( seq1part )
			lettercounts2 = Counter( seq2part )

			seq1gapless = seq1part.replace("-","")
			seq2gapless = seq2part.replace("-","")

			if not seq1gapless and not seq2gapless:
				merge_string += seq1part
				outcome_counter_dict["empty"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, both seqs are empty".format(pair, part, taxon1)
				continue

			# if parts are identical, take one and move on
			if seq1part == seq2part:
				merge_string += seq1part
				outcome_counter_dict["both-same"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, both seqs are identical".format(pair, part, taxon1)
				continue

			# parts are not the same, so if one seq is empty, always take the other
			if lettercounts1["-"] == partlength:
				merge_string += seq2part
				outcome_counter_dict["1st-none-2nd-yes"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, {} is empty".format(pair, part, taxon2, taxon1)
				continue
			elif lettercounts2["-"] == partlength: # even if the other is also empty
				merge_string += seq1part
				outcome_counter_dict["1st-yes-2nd-none"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, {} is empty".format(pair, part, taxon1, taxon2)
				continue

			### parts may or may not be overlapping
			# if one is full length and the other is not, take the full one
			if len(seq1gapless)==partlength and len(seq2gapless) < partlength:
				merge_string += seq1part
				outcome_counter_dict["1st-yes-2nd-partial"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, {} is full-length".format(pair, part, taxon1, taxon1)
				continue
			elif len(seq2gapless)==partlength and len(seq1gapless) < partlength:
				merge_string += seq2part
				outcome_counter_dict["1st-partial-2nd-yes"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {}, {} is full-length".format(pair, part, taxon2, taxon2)
				continue
			# if both are full length but different  #19066
			elif len(seq1gapless)==partlength and len(seq1gapless)==len(seq2gapless):
				merge_string += seq1part
				outcome_counter_dict["both-different"] += 1
				if verbose:
					print >> sys.stderr, "For {}: Part {}: using {} by default, full-length seqs are different".format(pair, part, taxon1)
				continue

			# if both are partial, merge the two sequences  #4491 #15944 #20055
			# in conflict, go with the letter of the longer seq
			# if seq length is equal, go with seq 1
			if len(seq1gapless)==len(seq2gapless) or len(seq1gapless) > len(seq2gapless):
				letterindex = 0
			else:
				letterindex = 1
			part_merge = ""
			for l1,l2 in izip(seq1part, seq2part):
				if l1 == "-": # l1 is gap, always take l2
					part_merge += l2
				elif l2 == "-": # l1 is not gap, l2 is gap
					part_merge += l1
				elif l1==l2: # letters are the same
					part_merge += l1
				else: # letters are neither gaps nor the same
					part_merge += [l1,l2][letterindex]
			merge_string += part_merge
			outcome_counter_dict["merged-partials"] += 1
			if verbose:
				print >> sys.stderr, "For {}: Part {}: using merged part, where {} was dominant".format(pair, part, [taxon1,taxon2][letterindex] )

		print >> sys.stdout, ">{}\n{}".format(merge_id, merge_string)
		for k in sorted( outcome_counter_dict.keys() ):
			print >> sys.stderr, "{}\t{}\t{}".format(pair, k, outcome_counter_dict[k] )

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-t','--taxa-pairs', nargs="*", help="space-separated list of pairs, divided by comma")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	parser.add_argument('-v','--verbose', action="store_true", help="additional output to log")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	merge_pairs(args.alignment, args.format, partitions, args.taxa_pairs, args.verbose)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
