#!/usr/bin/env python
#
# read_logs_for_hmm_stats.py created 2018-09-04
# v1.1 2023-01-20 python3 update

'''read_logs_for_hmm_stats.py  last modified 2023-01-20
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
	sys.stderr.write( "# reading gene names from {}\n".format(pairstatfile) )
	for line in open(pairstatfile,'r'):
		line = line.strip()
		lsplits = line.split()
		if lsplits[0]=="partition": # check for header line
			continue
		genehit = lsplits[1]
		gene_name = genehit.split("|")[-1].split("_")[0]
		genelist.append(gene_name)
	sys.stderr.write( "# found {} gene names\n".format( len(genelist) ) )
	return genelist

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-p', '--pair-stats', help="pair stats file with seq names in the same order")
	parser.add_argument('-f', '--fasta', help="fasta files of contigs or scaffolds")
	parser.add_argument('-l', '--log', nargs='*', help="log files from previous runs, in verbose mode", required=True)
	parser.add_argument('-n', '--names', nargs='*', help="names of previous runs, such as version name or letters")
	args = parser.parse_args(argv)

	categories = ["n_taxa", "evalue", "bpl"]
	seqcount_dict = defaultdict(list) # key is name, value is list of taxa counts
	evalue_dict = defaultdict(list) # key is name, value is list of evalues
	bpl_dict = defaultdict(list) # key is name, value is list of bits-per-length

	genelist = read_pair_stats(args.pair_stats)

	if len(args.log) != len(args.names):
		sys.exit("ERROR: {} log files for {} names".format( len(args.log), len(args.names) ) )

	for logfile, version_name in zip(args.log, args.names):
		linecounter = 0
		seqcounter = 0
		evaluecounter = 0
		bplcounter = 0
		sys.stderr.write( "# reading {} log information from {}\n".format( version_name, logfile ) )
		for line in open(logfile, 'r'):
			linecounter += 1
			if line[0]=="#": # information for hmm stats is in comments
				if line.find("# Found")==0: # n seqs in self-hmm
					ntaxa = line.strip().split(' ')[2]
					seqcounter += 1
					seqcount_dict[version_name].append(ntaxa)
				if line.find("# calculated e-value")==0:
					# calculated e-value for /mnt/data/est/supermatrix/v12/20180813-125832_partitions/v11_linsi_c10trim_no_breaks_6947_7796_part.aln as 1.000e-300
					evalue = line.strip().split(' ')[6]
					evaluecounter += 1
					evalue_dict[version_name].append(evalue)
				if line.find("# calculated bits")==0:
					# calculated bits-per-length cutoff for /mnt/data/est/supermatrix/v12/20180813-125832_partitions/v11_linsi_c10trim_no_breaks_6947_7796_part.aln as 1.623
					bpl = line.strip().split(' ')[7]
					bplcounter += 1
					bpl_dict[version_name].append(bpl)
				if line.find("# using minimum bits")==0:
					# using minimum bits-per-length cutoff for /mnt/data/est/supermatrix/v12/20180813-125832_partitions/v11_linsi_c10trim_no_breaks_7797_7960_part.aln of 0.90 Mon Aug 13 12:58:54 2018
					bpl = line.strip().split(' ')[8]
					bplcounter += 1
					bpl_dict[version_name].append(bpl)

		sys.stderr.write( "# found {} evalues, {} BPLs, with {} hmm-refs in {}\n".format( len(evalue_dict[version_name]), len(bpl_dict[version_name]), len(seqcount_dict[version_name]), version_name ) )

	name_cats = []
	for versionname in args.names:
		for cat in categories:
			name_cats.append("{}_{}".format(versionname, cat) )
	headerline = "gene_name\t{}\n".format( "\t".join(name_cats) )
	sys.stdout.write( headerline )
	for i, gene in enumerate(genelist):
		scorelist = []
		for version_name in args.names:
			scorelist.append( seqcount_dict[version_name][i] )
			scorelist.append( evalue_dict[version_name][i] )
			scorelist.append( bpl_dict[version_name][i] )
		geneline = "{}\t{}\n".format( gene, "\t".join(scorelist) )
		sys.stdout.write( geneline )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
