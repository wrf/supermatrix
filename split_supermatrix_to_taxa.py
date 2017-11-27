#!/usr/bin/env python
#
# split_supermatrix_to_taxa.py created 2017-09-18

'''split_supermatrix_to_taxa.py v1.0 2017-11-20
tool to split supermatrices into fasta files for each taxa

split_supermatrix_to_taxa.py -a matrix.phy -p partitions.txt -d taxa_dir -f phylip-relaxed

    for partitioned alignments, formats (-f) include:
  clustal, fasta, nexus, phylip, phylip-relaxed, stockholm

    large matrices can be gzipped, as .gz
'''

import sys
import os
import argparse
import time
import gzip
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

def split_taxa(fullalignment, alignformat, taxadir, partitions):
	'''read large alignment and write one fasta file for each taxa containing all proteins'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	print >> sys.stderr, "# alignment contains {} taxa with {} sites".format(len(alignedseqs), alignedseqs.get_alignment_length() ), time.asctime()
	filecounter = 0
	for seqrec in alignedseqs:
		taxafilename = "{}.fasta".format( os.path.join(taxadir, seqrec.id) )
		with open(taxafilename,'w') as tf:
			for part in partitions:
				seqpart = seqrec.seq[part[0]-1:part[1]] # correct alignment positions to python index
				if str(seqpart).replace("-","").replace("X",""): # ignore completely empty entries
					newseq = SeqRecord( seqrec.seq[part[0]-1:part[1]] , id="{}_{}-{}".format(seqrec.id, *part), description="" )
					tf.write( newseq.format("fasta") )
		filecounter += 1
	print >> sys.stderr, "# split alignment into {} files".format(filecounter), time.asctime()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-d','--taxa-directory', default="./", help="optional directory for taxa files")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	if not os.path.isdir(args.taxa_directory):
		if os.path.isfile(args.taxa_directory):
			raise OSError("ERROR: Cannot create directory {}, exiting".format(args.taxa_directory) )
		else:
			os.mkdir(args.taxa_directory)

	partitions = get_partitions(args.partition)
	split_taxa(args.alignment, args.format, args.taxa_directory, partitions)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
