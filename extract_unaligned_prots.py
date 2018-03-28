#!/usr/bin/env python
#
# extract_unaligned_prots.py  created by WRF 2018-03-21

'''extract_unaligned_prots.py  last modified 2018-03-22
  extract unaligned proteins from previous add_taxa_to_align runs and
  write each partition to a single file

extract_unaligned_prots.py -p partitions.txt -a */*

  -a should be a directory level where echo *part.fasta will find files
  for example, if files are nested as:

echo */*/*_part.fasta

  then -a should be */*/

  to then align all sequences with MAFFT:
for FILE in *.fasta; do BASE="${FILE%.fasta}"; mafft-linsi $FILE > $BASE.aln ; done
'''

import sys
import os
import argparse
import time
import glob
from collections import defaultdict
from Bio import SeqIO

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

def make_glob_list(part, alignmentfolder):
	'''from each partition, return the globbed list of filenames'''
	globstring = "{}*{}_{}_part.fasta".format(alignmentfolder, *part)
	#print >> sys.stderr, globstring
	globlist = glob.glob(globstring)
	return globlist

def extract_combine_fasta(fastalist, part, newdirectory):
	'''for each part, make a new fasta file containing all of the sequences for this part'''
	seqcounter = 0
	namecounter = defaultdict(int)
	combinedfasta = os.path.join(newdirectory, "unaligned_{}_{}_part.fasta".format(*part))
	print >> sys.stderr, "# generating file for partition {}".format(part), time.asctime()
	with open(combinedfasta,'w') as cf:
		for fastafile in fastalist:
			for seqrec in SeqIO.parse(fastafile, "fasta"):
				seqcounter += 1
				seqid = str(seqrec.id)
				if seqid in namecounter:
					namecounter[seqid] += 1
					seqid = "{}{}".format(seqid, namecounter[seqid])
					seqrec.id = seqid
					seqrec.description = ""
				else:
					namecounter[seqid] += 1
				cf.write(seqrec.format("fasta"))
	print >> sys.stderr, "# found {} sequences for partition {}".format(seqcounter, part), time.asctime()

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignments', help="folder tree structure for globbing")
	parser.add_argument('-d','--directory', default="./", help="optional directory for extracted files")
	parser.add_argument('-p','--partition', help="partition file for splitting large alignments")
	args = parser.parse_args(argv)

	if not os.path.isdir(args.directory):
		if os.path.isfile(args.directory):
			raise OSError("ERROR: Cannot create directory {}, exiting".format(args.directory) )
		else:
			os.mkdir(args.directory)

	partitions = get_partitions(args.partition)

	for part in partitions:
		fastalist = make_glob_list(part, args.alignments)
	#	print >> sys.stdout, fastalist
		extract_combine_fasta(fastalist, part, args.directory)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
