#!/usr/bin/env python
#
# reorder_matrix_by_cov.py created 2017-10-24

'''reorder_matrix_by_cov.py  last modified 2018-02-27

reorder_matrix_by_cov.py -a matrix.phy -p partitions.txt -o reordered_matrix.phy

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    to take the highest coverage genes stopping after N sites, use -m
reorder_matrix_by_cov.py -a matrix.phy -p partitions.txt -o reordered_matrix.phy -m 5000

    add gene names in RAxML-type partition format, using --pair-stats
    pair stats file is the output of align_pair_stats.py, or any file as:
1-1000    GENE
    or can include species or SwissProt ID
Genus_species_1-1000    sp|123456|GENE
'''

import sys
import os
import argparse
import time
import gzip
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

def parts_to_genes(pairstatsfile):
	'''read tabular pair-wise gene stats, return a dict where key is partition and value is gene'''
	part_to_gene = {}
	print >> sys.stderr, "# reading partitions and gene names from {}".format(pairstatsfile), time.asctime()
	for line in open(pairstatsfile,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			if lsplits[0]=="partition": # skip first line
				continue
			gene = lsplits[1].split("|")[-1]
			partition = lsplits[0].split("_")[-1]
			part_to_gene[partition] = gene
	print >> sys.stderr, "# found {} gene names".format(len(part_to_gene)), time.asctime()
	return part_to_gene

def reorder_alignments(fullalignment, alignformat, partitions, maxsites, genenamedict):
	'''read alignment, and return a new alignment where partitions are reordered from highest coverage to lowest'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()

	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	num_species = len(alignedseqs)
	newalign = alignedseqs[:,0:0] # start with blank alignment
	newpartitions = []
	print >> sys.stderr, "# alignment has {} taxa with {} positions".format(num_species, alignedseqs.get_alignment_length()), time.asctime()

	newnamedict = {} # key is new partition, value is gene name from old partition
	alignparts = {} # key is old partition, value is alignment part
	aligncompleteness = {} # key is old partition, value is coverage
	print >> sys.stderr, "# scoring partitions", time.asctime()
	for part in partitions:
		alignpart = alignedseqs[:, part[0]-1:part[1] ] # alignment of each partition only
		gaplist = []
		#completecounter = {0:0, 1:0, 2:0}
		for i,seqrec in enumerate(alignpart):
			species = seqrec.id
			seqlen = len(seqrec.seq)
			lettercounts = Counter(str(seqrec.seq).replace("X","-"))
			gapcount = lettercounts.get("-",0) + lettercounts.get("?",0)
			occupancyscore = 100 - 100*gapcount/seqlen
			gaplist.append(occupancyscore)
			#occupancyscore = 2 # by default is present, reassign if absent or partial
			#if lettercounts["-"] == seqlen or lettercounts["?"] == seqlen: # seq is all gaps, so no seq
			#	occupancyscore = 0 # set to 0 if all gaps
			#elif lettercounts["-"] >= seqlen * 0.8: # partial means 20% or more of sequence is gaps
			#	occupancyscore = 1 # set to 1 if partial
			#completecounter[occupancyscore] += 1
		#covscore = sum( k*v for k,v in completecounter.items() )
		covscore = sum( gaplist )
		#print >> sys.stderr, part, covscore
		aligncompleteness[part] = covscore
		alignparts[part] = alignpart

	print >> sys.stderr, "# sorting partitions", time.asctime()
	for part, score in sorted(aligncompleteness.items(), key=lambda x: x[1], reverse=True):
		#print >> sys.stderr, part, score
		newindices = "{}:{}".format( newalign.get_alignment_length()+1, newalign.get_alignment_length()+part[1]-part[0]+1 )
		newpartitions.append( newindices )
		newalign += alignparts[part]
		if genenamedict: # transfer name to new indices
			newnamedict[newindices] = genenamedict.get("{}-{}".format(*part),None)
		if maxsites is not None and newalign.get_alignment_length() >= maxsites:
			break
	print >> sys.stderr, "# alignment has {} partitions with {} sites".format( len(newpartitions), newalign.get_alignment_length() ), time.asctime()
	return newalign, newpartitions, newnamedict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-m','--max-sites', type=int, help="optional limit to stop after N sites")
	parser.add_argument('-o','--output', help="output name for new alignment", required=True)
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	parser.add_argument('--pair-stats', help="pair stats file for gene names, to use with RAxML style partition format")
	args = parser.parse_args(argv)

	partitions = get_partitions(args.partition)
	genenames = parts_to_genes(args.pair_stats) if args.pair_stats else None
	sortedalignment, partitionlist, genenames = reorder_alignments(args.alignment, args.format, partitions, args.max_sites, genenames)
	AlignIO.write(sortedalignment, args.output, args.format)
	print >> sys.stderr, "# Supermatrix written to {}".format(args.output), time.asctime()
	with open("{}.partition.txt".format(args.output),'w') as pf:
		if args.pair_stats:
			for part in partitionlist:
				print >> pf, "model, {} = {}-{}".format( genenames[part], *part.split(":") )
		else:
			print >> pf, ",".join(partitionlist)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
