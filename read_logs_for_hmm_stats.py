#!/usr/bin/env python
#
# read_logs_for_hmm_stats.py created 2018-09-04

'''read_logs_for_hmm_stats.py  last modified 2018-09-04
    extract hmm stats for each log file in a series of runs

read_logs_for_hmm_stats.py -p pair_stats.tab -l v1.log v2.log -n v1 v2 > hmm_stats.tab

    pair stats file is tab-delimited, where first column is partitions,
    second column is gene name, or SwissProt name
partition	protID	trimmedLength	trimmedPercent	span	spanLength	spanPercent	refProtLength
Homo_sapiens_1-136	sp|O15145|ARPC3_HUMAN	136	0.764	(28, 177)	149	0.837	178
'''

import sys
import argparse
from collections import defaultdict

def read_pair_stats(pairstatfile):
	'''read pair stats and return a list of gene names'''
	genelist = []
	print >> sys.stderr, "# reading gene names from {}".format(pairstatfile)
	for line in open(pairstatfile,'r'):
		line = line.strip()
		lsplits = line.split()
		if lsplits[0]=="partition": # check for header line
			continue
		genehit = lsplits[1]
		gene_name = genehit.split("|")[-1].split("_")[0]
		genelist.append(gene_name)
	print >> sys.stderr, "# found {} gene names".format( len(genelist) )
	return genelist

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-p', '--pair-stats', help="pair stats file with seq names in the same order")
	parser.add_argument('-f', '--fasta', help="fasta files of contigs or scaffolds")
	parser.add_argument('-l', '--log', nargs='*', help="log files from previous runs, in verbose mode")
	parser.add_argument('-n', '--names', nargs='*', help="names of previous runs, such as version name or letters")
	args = parser.parse_args(argv)

	evalue_dict = defaultdict(list) # key is name, value is list of evalues

	genelist = read_pair_stats(args.pair_stats)

	if len(args.log) != len(args.names):
		sys.exit("ERROR: {} log files for {} names".format( len(args.log), len(args.names) ) )

	for logfile, version_name in zip(args.log, args.names):
		linecounter = 0
		evaluecounter = 0
		print >> sys.stderr, "# reading {} log information from {}".format( version_name, logfile )
		for line in open(logfile, 'r'):
			linecounter += 1
			if line[0]=="#": # information for hmm stats is in comments
				if line.find("# calculated e-value")==0:
					evalue = line.strip().split(' ')[6]
					evaluecounter += 1
					evalue_dict[version_name].append(evalue)
		print >> sys.stderr, "# found evalues for {} genes in {}".format( len(evalue_dict[version_name]), version_name )
	#	print >> sys.stderr, evalue_dict

	headerline = "gene_name\t{}".format( "\t".join(args.names) )
	print >> sys.stdout, headerline
	for i, gene in enumerate(genelist):
		evaluelist = []
		for version_name in args.names:
			evaluelist.append( evalue_dict[version_name][i] )
		geneline = "{}\t{}".format( gene, "\t".join(evaluelist) )
		print >> sys.stdout, geneline

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
