#!/usr/bin/env python
#
# split_supermatrix_to_genes.py created 2018-02-12

'''split_supermatrix_to_genes.py v1.0 2018-02-14
tool to re-partition supermatrices into fasta files for each gene

split_supermatrix_to_genes.py -a matrix.phy -p partitions.txt -d aln_dir -f phylip-relaxed

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

def split_genes(fullalignment, alignformat, alndir, partitions, fileprefix):
	'''read large alignment and write one fasta file for each taxa containing all proteins'''
	if fullalignment.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading alignment {} as gzipped".format(fullalignment), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading alignment {}".format(fullalignment), time.asctime()
	alignedseqs = AlignIO.read(opentype(fullalignment), alignformat)
	num_species = len(alignedseqs)
	al_length = alignedseqs.get_alignment_length()
	print >> sys.stderr, "# alignment contains {} taxa with {} sites".format( num_species, al_length ), time.asctime()
	filecounter = 0
	for part in partitions:
		alignpart = alignedseqs[:, part[0]-1:part[1] ] # alignment of each partition only
		maxdigits = str(len(str(al_length)))
		bufferedpartstring = fileprefix + "_{:0"+maxdigits+"}_{:0"+maxdigits+"}.aln"
		genefilename = "{}".format( os.path.join(alndir, bufferedpartstring.format(*part)) )
		AlignIO.write(alignpart, genefilename, "fasta")
		filecounter += 1
	print >> sys.stderr, "# split alignment into {} files".format(filecounter), time.asctime()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="supermatrix alignment")
	parser.add_argument('-d','--gene-directory', default="./", help="optional directory for alignment files")
	parser.add_argument('-f','--format', default="fasta", help="alignment format [fasta]")
	parser.add_argument('-g','--gene-prefix', default="gene", help="prefix for each file name [gene]")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	if not os.path.isdir(args.gene_directory):
		if os.path.isfile(args.gene_directory):
			raise OSError("ERROR: Cannot create directory {}, exiting".format(args.gene_directory) )
		else:
			os.mkdir(args.gene_directory)

	partitions = get_partitions(args.partition)
	split_genes(args.alignment, args.format, args.gene_directory, partitions, args.gene_prefix)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
